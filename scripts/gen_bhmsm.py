USE_VERSION = 'v0.2b'
import os
import sys
import yt
import caesar
import unyt
import argparse
import numpy as np
import pandas as pd
from copy import copy, deepcopy
from scipy.optimize import curve_fit

import gen_sim_data as gen

from velociraptor.observations.objects import ObservationalData, MultiRedshiftObservationalData
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

##########################################################################################################################
#### This function generates galaxy stellar mass functions (GSMF) from caesar files for the given simulation in a given directory.
#### Call it using python gen_bhmsm.py --modeldir=path/to/model --simdir=SIMDIR --sim=SIM
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



# Set units and labels
x_units = unyt.Solar_Mass
y_units = unyt.Solar_Mass

x_label = r'$\log(M_{\star}\,/\,%s)$' % x_units.units.latex_representation()
y_label = r'$\log(M_{BH}\,/\,%s)$' % x_units.units.latex_representation()


# Min and max stellar masses
min_log_gal_stellar_mass = 8
max_log_gal_stellar_mass = 14

min_log_gal_stellar_mass_for_fitting = 10
max_log_gal_stellar_mass_for_fitting = 13


## Minimum number of galaxies required in a bin when binning
minN = 10

## x spacing for binning
dx = 0.25




def gen_observable(SNAPLIST, SNAPDIR, SIM):
    print (SNAPDIR)
    '''
    SNAPLIST: list of the snapshot numbers
    SIM: Sim name, e.g. m50n512, used in snapshot name (assumes Simba naming convention)
    SNAPDIR = Directory where snapshots are located
    '''
    for j in SNAPLIST:
        print(j)
        # path to the file
        SNAP  = os.path.join(SNAPDIR, '%s_%04d.hdf5' % (SIM, j))
        if not os.path.exists(SNAP):
            print(SNAP, "does not exist")
            continue

        CAESARFILE = os.path.join(SNAPDIR, '%s_caesar_%04d.hdf5' % (SIM, j))

        ds = yt.load(SNAP)
        Vcom = ds.domain_width.in_units("Mpccm").prod()  # Comoving box volume in Mpc^3

        print('Loading caesar object')
        obj = caesar.load(CAESARFILE)
        redshift = obj.simulation.redshift
        h0 = obj.simulation.hubble_constant

        print("Comoving Volume=%s" % Vcom)
        print('z=%s' % (redshift))
        print('h0=%s' % (h0))


        # Cosmological Parameters
        omega_matter = ds.cosmology.omega_matter
        cosmology = FlatLambdaCDM(name=r'Flat $\Lambda$CDM', H0=h0*100. * u.km / u.s / u.Mpc, Om0=omega_matter)


        
        # Dictionary to save all observables in
        save_data = {}

        

        ## Simulation data for observable
        gal_stellar_masses = np.array([gal.masses['stellar'].in_units(x_units) for gal in obj.galaxies])
        gal_bh_masses = np.array([gal.masses['bh'].in_units(x_units) for gal in obj.galaxies])
        
        log_gal_stellar_masses = np.log10(gal_stellar_masses)
        log_gal_bh_masses = np.log10(gal_bh_masses)
        

        gal_stellar_masses_no_zeros = gal_stellar_masses[gal_bh_masses > 0]
        gal_bh_masses_no_zeros = gal_bh_masses[gal_bh_masses > 0]
        
        log_gal_stellar_masses_no_zeros = np.log10(gal_stellar_masses_no_zeros)
        log_gal_bh_masses_no_zeros = np.log10(gal_bh_masses_no_zeros)
        
        
        gal_stellar_masses_for_fitting = gal_stellar_masses_no_zeros[
            np.where(np.logical_and(
                gal_stellar_masses_no_zeros >= 10**min_log_gal_stellar_mass_for_fitting, 
                gal_stellar_masses_no_zeros <= 10**max_log_gal_stellar_mass_for_fitting))]
        gal_bh_masses_for_fitting = gal_bh_masses_no_zeros[
            np.where(np.logical_and(
                gal_stellar_masses_no_zeros >= 10**min_log_gal_stellar_mass_for_fitting, 
                gal_stellar_masses_no_zeros <= 10**max_log_gal_stellar_mass_for_fitting))]
        
        log_gal_stellar_masses_for_fitting = np.log10(gal_stellar_masses_for_fitting)
        log_gal_bh_masses_for_fitting = np.log10(gal_bh_masses_for_fitting)


        save_data['log_data'] = {
            'x':log_gal_stellar_masses_no_zeros * x_units,
            'xerr':np.zeros(len(log_gal_stellar_masses_no_zeros)) * x_units,
            'y':log_gal_bh_masses_no_zeros * y_units,
            'yerr':np.zeros(len(log_gal_bh_masses_no_zeros)) * y_units, 
            'x_label':x_label, 
            'y_label':y_label,
            'plot_as':'points',
        }
    
    
        
        ## Bin data and find means/medians in each bin
        x_bin_centres, y_bin_means, y_bin_stds, y_bin_medians, y_bin_quantiles_lo, y_bin_quantiles_hi = gen.regular_bin(
            log_gal_stellar_masses_no_zeros, log_gal_bh_masses_no_zeros, 
            dx=dx, min_x=min_log_gal_stellar_mass, max_x=max_log_gal_stellar_mass, 
            calc_min_x=False, calc_max_x=False, minN=minN)
        
        save_data['binned_log_mean'] = {
            'x':x_bin_centres * x_units,
            'xerr':np.zeros(len(x_bin_centres)) * x_units,
            'y':y_bin_means * y_units,
            'yerr':y_bin_stds * y_units, 
            'x_label':x_label, 
            'y_label':y_label,
            'plot_as':'line',
        }
        save_data['binned_log_median'] = {
            'x':x_bin_centres * x_units,
            'xerr':np.zeros(len(x_bin_centres)) * x_units,
            'y':y_bin_medians * y_units,
            'yerr':(np.abs(y_bin_medians - y_bin_quantiles_lo), np.abs(y_bin_medians - y_bin_quantiles_hi)) * y_units, 
            'x_label':x_label, 
            'y_label':y_label,
            'plot_as':'line',
        }
        
        
        
        ## Variable binning such that there are EXACTLY X objects (cannot be set) in N bins (must be set)
#         x_bin_centres, y_bin_means, y_bin_stds, y_bin_medians, y_bin_quantiles_lo, y_bin_quantiles_hi = gen.exact_count_bin(
#             log_gal_stellar_masses_no_zeros, log_gal_bh_masses_no_zeros,
#             Nbins=20, min_x=min_log_gal_stellar_mass, max_x=max_log_gal_stellar_mass, 
#             calc_min_x=False, calc_max_x=False)
        
#         save_data['exact_count_binned_log_mean'] = {
#             'x':x_bin_centres * x_units,
#             'xerr':np.zeros(len(x_bin_centres)) * x_units,
#             'y':y_bin_means * y_units,
#             'yerr':y_bin_stds * y_units, 
#             'x_label':x_label, 
#             'y_label':y_label,
#             'plot_as':'line',
#         }
#         save_data['exact_count_binned_log_median'] = {
#             'x':x_bin_centres * x_units,
#             'xerr':np.zeros(len(x_bin_centres)) * x_units,
#             'y':y_bin_medians * y_units,
#             'yerr':(np.abs(y_bin_medians - y_bin_quantiles_lo), np.abs(y_bin_medians - y_bin_quantiles_hi)) * y_units, 
#             'x_label':x_label, 
#             'y_label':y_label,
#             'plot_as':'line',
#         }
        
        
        
        
        ## Variable binning such that there are AT LEAST N objects in each bin
#         x_bin_centres, y_bin_means, y_bin_stds, y_bin_medians, y_bin_quantiles_lo, y_bin_quantiles_hi = gen.min_count_bin(
#             log_gal_stellar_masses_no_zeros, log_gal_bh_masses_no_zeros,
#             x_lo=min_log_gal_stellar_mass, x_hi=max_log_gal_stellar_mass, 
#             calc_min_x=True, calc_max_x=False,
#             min_N_obj=3, min_bin_width=0.5, delta_x=0.02)
        
#         save_data['min_count_binned_log_mean'] = {
#             'x':x_bin_centres * x_units,
#             'xerr':np.zeros(len(x_bin_centres)) * x_units,
#             'y':y_bin_means * y_units,
#             'yerr':y_bin_stds * y_units, 
#             'x_label':x_label, 
#             'y_label':y_label,
#             'plot_as':'line',
#         }
#         save_data['min_count_binned_log_median'] = {
#             'x':x_bin_centres * x_units,
#             'xerr':np.zeros(len(x_bin_centres)) * x_units,
#             'y':y_bin_medians * y_units,
#             'yerr':(np.abs(y_bin_medians - y_bin_quantiles_lo), np.abs(y_bin_medians - y_bin_quantiles_hi)) * y_units, 
#             'x_label':x_label, 
#             'y_label':y_label,
#             'plot_as':'line',
#         }
        
        
        
        
        ## Fit from KH13 (pg. 57, eqn 10)
        ## They use a symmetric, least-squares ﬁt
        log_x_arr = np.linspace(min_log_gal_stellar_mass, max_log_gal_stellar_mass, 10)

        log_x_arr, log_y, log_y_err = gen.log_fit_fn_with_errors_for_fitting(
            log_x_arr, log_gal_stellar_masses_for_fitting, log_gal_bh_masses_for_fitting, 
            func=None, paper='kh13', determine_func=True)
        
        save_data['kh13_log_fit'] = {
            'x':log_x_arr * x_units,
            'xerr':np.zeros(len(log_x_arr)) * x_units,
            'y':log_y * y_units,
            'yerr':log_y_err * y_units, 
            'x_label':x_label, 
            'y_label':y_label,
            'plot_as':'line',
        }

        
        


        # Save to HDF5 files
        output_directory = SNAPDIR

        if not os.path.exists(output_directory):
            os.mkdir(output_directory)

        for key, val in save_data.items():
            output_filename = "%s_bhmsm_%s_min%sgal_%04d.hdf5" % (SIM, key, minN, j)

            comment = f"h-corrected for SWIFT using cosmology: {cosmology}."
            citation = ""
            bibcode = ""
            name = r"$M_{BH}-M_*$ Relation"
            plot_as = val['plot_as']
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
    print('DONE')


                                              
gen_observable(SNAPLIST, os.path.join(MODELDIR, SIMDIR), SIM)

#print()
print("DONE")
print()
print()