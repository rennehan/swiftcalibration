USE_VERSION = 'v0.2b'
import os
import sys
import yt
import caesar
import unyt
import argparse
import numpy as np

import gen_sim_data as gen

from velociraptor.observations.objects import ObservationalData, MultiRedshiftObservationalData
#from astropy.cosmology import WMAP7 as cosmology
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

##########################################################################################################################
#### This function generates galaxy stellar mass functions (GSMF) from caesar files for the given simulation in a given directory.
#### Call it using python gen_gsmf.py --modeldir=path/to/model --simdir=SIMDIR --sim=SIM
#### Based on Renier Hough's "Create_caesarfiles.py"
##########################################################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('--modeldir')
parser.add_argument('--simdir')
parser.add_argument('--sim', type=str)
args = parser.parse_args()

MODELDIR = args.modeldir
SIMDIR = args.simdir
SIM = args.sim


## Need to set this in here
#SNAPLIST = list(range(0,16))
SNAPLIST = [15]

print (SNAPLIST)


## Units and labels
x_units = unyt.Solar_Mass
length_units = 'Mpccm'
y_units = unyt.Mpc**(-3)

x_label = r'$\log(M_{\star}\,/\,%s)$' % x_units.units.latex_representation()
y_label = r'$\log(\Phi=dn/dlogM_\star\,[\,\mathrm{dex}^{-1}\,\mathrm{cMpc}^{-3}])$'



# Set min and max stellar masses for binning
min_log_gal_stellar_mass = 6
max_log_gal_stellar_mass = 14


## Minimum number of galaxies require in a bin when binning
minN = 10

## x spacing for binning
dx = 0.2


def gen_observable(SNAPLIST, SNAPDIR, SIM):
    print (SNAPDIR)
    '''
    SNAPLIST: list of the snapshot numbers
    SIM: Sim name, e.g. m50n512, used in snapshot name (assumes Simba naming convention)
    SNAPDIR = Directory where snapshots are located
    '''
    for j in SNAPLIST:
        print (j)
        # path to the file
        SNAP  = os.path.join(SNAPDIR, '%s_%04d.hdf5' % (SIM,j))
        if not os.path.exists(SNAP):
            print(SNAP, "does not exist")
            continue

        CAESARFILE = os.path.join(SNAPDIR, '%s_caesar_%04d.hdf5' % (SIM,j))

        ds = yt.load(SNAP)
        Vcom = ds.domain_width.in_units(length_units).prod()  # Comoving box volume in Mpc^3

        print('Loading caesar object')
        obj = caesar.load(CAESARFILE)
        redshift = obj.simulation.redshift
        h0 = obj.simulation.hubble_constant
        boxsize = obj.simulation.boxsize.in_units(length_units)

        print("Comoving Volume=%s" % Vcom)
        print('z=%s' % (redshift))
        print('h0=%s' % (h0))


        # Cosmological Parameters
        omega_matter = ds.cosmology.omega_matter
        cosmology = FlatLambdaCDM(name=r'Flat $\Lambda$CDM', H0=h0*100. * u.km / u.s / u.Mpc, Om0=omega_matter)
        
        
        
        # Dictionary to save all observables in
        save_data = {}
        

            
        ## Stellar masses from simulation
        stellar_masses = np.array([gal.masses['stellar'].in_units(x_units) for gal in obj.galaxies])
        log_stellar_masses = np.log10(stellar_masses)
        pos_Mpc = np.array([gal.pos.in_units(length_units) for gal in obj.galaxies])
        
        
        # Add eddington bias to calculated/measured/observed stellar masses (for GSMF only?)
        # See Kugel+23 (FLAMINGO: Calibrating large cosmological hydrodynamical simulations with machine learning)
        # Section 3.1.1 eqn 11
        rng = np.random.default_rng()
        mu_eddington = 0
        sigma_eddington = min(0.070+0.071*redshift, 0.3) # dex
        # eddington_bias = rng.lognormal(mean=mu_eddington, sigma=sigma_eddington)#, size=1)
        eddington_bias = rng.normal(loc=mu_eddington, scale=sigma_eddington, size=len(log_stellar_masses))

        log_stellar_masses_eddington_biased = log_stellar_masses + eddington_bias
        stellar_masses_eddington_biased = 10**log_stellar_masses_eddington_biased
        
        
        
        
        ## Unbiased stellar mass function
        
        logM_ax, logPhi, logPhi_total_lo_err, logPhi_total_hi_err, logPhi_total_err, logPhi_total_err_v2 = gen.mass_function_with_error(
            stellar_masses, pos_Mpc, boxsize, Vcom, 
            dlogM=dx, min_logM=min_log_gal_stellar_mass, max_logM=max_log_gal_stellar_mass, 
            calc_min_logM=False, calc_max_logM=False, minN=minN)

        isfinite = np.isfinite(logPhi)
        logM_ax = logM_ax[isfinite]
        logPhi = logPhi[isfinite]
        logPhi_total_lo_err = logPhi_total_lo_err[isfinite]
        logPhi_total_hi_err = logPhi_total_hi_err[isfinite]
        logPhi_total_err = logPhi_total_err[isfinite]
        logPhi_total_err_v2 = logPhi_total_err_v2[isfinite]
    
    
        save_data['log_data'] = {
            'x':logM_ax * x_units,
            'xerr':np.zeros(len(logM_ax)) * x_units,
            'y':logPhi * y_units,
            'yerr':logPhi_total_err * y_units, 
            'x_label':x_label, 
            'y_label':y_label,
            'name':'GSMF',
        }
        
        
    
        
        ## Eddington biased stellar mass function
        
        logM_ax, logPhi, logPhi_total_lo_err, logPhi_total_hi_err, logPhi_total_err, logPhi_total_err_v2 = gen.mass_function_with_error(
            stellar_masses_eddington_biased, pos_Mpc, boxsize, Vcom, 
            dlogM=dx, min_logM=min_log_gal_stellar_mass, max_logM=max_log_gal_stellar_mass, 
            calc_min_logM=False, calc_max_logM=False, minN=minN)

        isfinite = np.isfinite(logPhi)
        logM_ax = logM_ax[isfinite]
        logPhi = logPhi[isfinite]
        logPhi_total_lo_err = logPhi_total_lo_err[isfinite]
        logPhi_total_hi_err = logPhi_total_hi_err[isfinite]
        logPhi_total_err = logPhi_total_err[isfinite]
        logPhi_total_err_v2 = logPhi_total_err_v2[isfinite]
    
    
        save_data['log_data_eddington_biased'] = {
            'x':logM_ax * x_units,
            'xerr':np.zeros(len(logM_ax)) * x_units,
            'y':logPhi * y_units,
            'yerr':logPhi_total_err * y_units, 
            'x_label':x_label, 
            'y_label':y_label,
            'name':'GSMF (Eddington Biased)',
        }
        
        
        
        
        # Save to Velociraptor HDF5 files
        output_directory = SNAPDIR

        if not os.path.exists(output_directory):
            os.mkdir(output_directory)

        for key, val in save_data.items():
            output_filename = "%s_gsmf_%s_min%sgal_%04d.hdf5" % (SIM, key, minN, j)

            comment = f"h-corrected for SWIFT using cosmology: {cosmology}."
            citation = ""
            bibcode = ""
            name = val['name']
            plot_as = "line"
            redshift_lower = redshift
            redshift_upper = redshift


            x = val['x']
            xerr = val['xerr']

            y = val['y']
            yerr = val['yerr']

            processed = ObservationalData()
            processed.associate_x(x, scatter=xerr, comoving=True, description=val['x_label'])
            processed.associate_y(y, scatter=yerr, comoving=True, description=val['y_label'])
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
            
            
        print()
        
    print()
    print()




gen_observable(SNAPLIST, os.path.join(MODELDIR, SIMDIR), SIM)

#print()
print("DONE")
print()
print()