USE_VERSION = 'v0.2b'
import os
import sys
import yt
import caesar
import argparse
import numpy as np

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

def gen_gsmf(SNAPLIST, SNAPDIR, SIM='doug_s18n128'):
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
    h = obj.simulation.hubble_constant
#    Vcom *= h**(-3)

    print("Comoving Volume=%s" % Vcom)
    print('z=%s' % (redshift))
    print('h=%s' % (h))



    stellar_masses = np.array([gal.masses['stellar'] for gal in obj.galaxies])

    logM  = np.log10(stellar_masses)        #Take logarithm
    dlogM = 0.5                             #bin width in dex
#    bins = np.arange(logM.min(), logM.max(), dlogM)
    bins = np.arange(8.0, 13.5, dlogM)
    nbins = len(bins)
#    nbins = 10                              #Number of bins to divide data into
    Phi,edg = np.histogram(logM,bins=bins) #Unnormalized histogram and bin edges
#    dlogM    = edg[1] - edg[0]                 #Bin size
    logM_ax   = edg[0:-1] + dlogM/2.               #Mass axis
    Phi   = Phi / Vcom / dlogM                 #Normalize to volume and bin size
    
#    print(Max)
#    print(Phi)
#    print(np.column_stack((Max, Phi)))
    np.savetxt(GSMFFILE, np.column_stack((logM_ax, Phi)))




for i in range(0, len(SIM)):
#    gen_gsmf(SNAPLIST, SNAPDIR, NBINS=NBINS, SIM=SIM[i])    
    gen_gsmf(SNAPLIST, os.path.join(SNAPDIR, SIM[i]))
