USE_VERSION = 'v0.2b'
import os
import sys
import yt
import caesar
import unyt
import argparse
import numpy as np

from velociraptor.observations.objects import ObservationalData, MultiRedshiftObservationalData
#from astropy.cosmology import WMAP7 as cosmology
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

##########################################################################################################################
#### This function generates galaxy stellar mass functions (GSMF) from ceasar files for the given simulation in a given directory.
#### Call it using python gen_gsmf.py SNAPDIR MODEL1 MODEL2...
#### Based on Renier Hough's "Create_caesarfiles.py"
##########################################################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('myArgs', nargs='*')
args = parser.parse_args()
SNAPDIR = args.myArgs[0]
#NBINS = int(args.myArgs[1])
SIM = []
for i in range(1,len(args.myArgs)):
   SIM.append(args.myArgs[i])

#SNAPLIST = list(range(0,16))
SNAPLIST = [15]

print (SNAPLIST)

def gen_gsmf(SNAPLIST, SNAPDIR, SIM='simba_s12.5n128'):
  print (SNAPDIR)
  '''
  SNAPLIST: list of the snapshot numbers
  SIM: Sim name, e.g. m50n512, used in snapshot name (assumes Simba naming convention)
  SNAPDIR = Directory where snapshots are located
  NBINS = number of bins with which to create GSMF
  CAESARLOC = Subdirectory to write Caesar catalogs
  '''
  for j in SNAPLIST:
    print (j)
    # path to the file
    SNAP  = os.path.join(SNAPDIR, '%s_%04d.hdf5' % (SIM,j))
    if not os.path.exists(SNAP):
      print(SNAP, "does not exist")
      continue
    
    CAESARFILE = os.path.join(SNAPDIR, '%s_caesar_%04d.hdf5' % (SIM,j))
    GSMFFILE = os.path.join(SNAPDIR, '%s_gsmf_%04d.txt' % (SIM,j))

    ds = yt.load(SNAP)
    Vcom = ds.domain_width.in_units("Mpccm").prod()  # Comoving box volume in Mpc^3
    
    print('Loading caesar object')
    obj = caesar.load(CAESARFILE)
    redshift = obj.simulation.redshift
    h0 = obj.simulation.hubble_constant
#    Vcom *= h**(-3)

    print("Comoving Volume=%s" % Vcom)
    print('z=%s' % (redshift))
    print('h0=%s' % (h0))
    
    
    # Cosmological Parameters
    omega_matter = ds.cosmology.omega_matter
    cosmology = FlatLambdaCDM(name=r'Flat $\Lambda$CDM', H0=h0*100. * u.km / u.s / u.Mpc, Om0=omega_matter)


    stellar_masses = np.array([gal.masses['stellar'] for gal in obj.galaxies])

    logM  = np.log10(stellar_masses)        #Take logarithm
    dlogM = 0.5                             #bin width in dex
#    bins = np.arange(logM.min(), logM.max(), dlogM)
    bins = np.arange(8.0, 13.5+dlogM, dlogM)
    nbins = len(bins)
#    nbins = 10                              #Number of bins to divide data into
    Phi, edg = np.histogram(logM,bins=bins) #Unnormalized histogram and bin edges
    Phi_err = np.sqrt(Phi)

#    dlogM    = edg[1] - edg[0]                 #Bin size
    logM_ax   = edg[0:-1] + dlogM/2.               #Mass axis

    Phi   = Phi / Vcom / dlogM                 #Normalize to volume and bin size
    Phi_err /= Vcom * dlogM
    
#    print(Max)
#    print(Phi)
#    print(np.column_stack((Max, Phi)))


    # Save to text file
    if os.path.exists(GSMFFILE):
       os.remove(GSMFFILE)
    
    np.savetxt(GSMFFILE, np.column_stack((10**logM_ax, Phi, Phi_err)), header='M/Msun\tPhi/cMpc^-3\tPhi_err/cMpc^-3')
    
    
    
    # Create and save Velociraptor HDF5 data file
    output_filename = "%s_gsmf_%04d.hdf5" % (SIM,j)
    output_directory = SNAPDIR

    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    processed = ObservationalData()

    comment = f"h-corrected for SWIFT using cosmology: {cosmology.name}."
    citation = ""
    bibcode = ""
    name = "GSMF from swimba"
    plot_as = "points"
#    redshift = redshift
    redshift_lower = redshift
    redshift_upper = redshift
    h = cosmology.h
    #print(h)

    log_M = logM_ax
    #M = 10 ** (log_M) * unyt.Solar_Mass / h
    M = 10 ** (log_M) * unyt.Solar_Mass
    #Phi = (10**raw.T[1] * (h ** 3)) * unyt.Mpc ** (-3)
    Phi = Phi.value * unyt.Mpc**(-3) #* unyt.mag ** (-1)

    Phi_err = Phi_err.value * unyt.Mpc**(-3)
    #print(Phi_err)

    processed.associate_x(M, scatter=None, comoving=True, description="Galaxy Stellar Mass")
    processed.associate_y(Phi, scatter=Phi_err, comoving=True, description="Phi (GSMF)")
    processed.associate_citation(citation, bibcode)
    processed.associate_name(name)
    processed.associate_comment(comment)
    processed.associate_redshift(redshift, redshift_lower, redshift_upper)
    processed.associate_plot_as(plot_as)
    processed.associate_cosmology(cosmology)

    output_path = os.path.join(output_directory, output_filename)

    if os.path.exists(output_path):
       os.remove(output_path)

    processed.write(filename=output_path)




for i in range(0, len(SIM)):
#    gen_gsmf(SNAPLIST, SNAPDIR, SIM=SIM[i])    
    gen_gsmf(SNAPLIST, os.path.join(SNAPDIR, SIM[i]))
