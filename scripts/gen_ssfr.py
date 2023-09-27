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
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

##########################################################################################################################
#### This function generates sSFR histogram from caesar files for the given simulation in a given directory.
#### Call it using python gen_ssfr_hist.py --modeldir=path/to/model --simdir=SIMDIR --sim=SIM
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

#SNAPLIST = list(range(0,16))
SNAPLIST = [15]

print (SNAPLIST)



mass_units = unyt.Solar_Mass
ssfr_units = unyt.Gyr**(-1)



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
        
        # Path to the snap's caesarfile
        CAESARFILE = os.path.join(SNAPDIR, '%s_caesar_%04d.hdf5' % (SIM,j))

        ds = yt.load(SNAP)
        Vcom = ds.domain_width.in_units("Mpccm").prod()  # Comoving box volume in Mpc^3

        print('Loading caesar object')
        obj = caesar.load(CAESARFILE)
        redshift = obj.simulation.redshift
        h0 = obj.simulation.hubble_constant
        boxsize = obj.simulation.boxsize.in_units("Mpccm")

        print("Comoving Volume=%s" % Vcom)
        print('z=%s' % (redshift))
        print('h0=%s' % (h0))


        # Cosmological Parameters
        omega_matter = ds.cosmology.omega_matter
        cosmology = FlatLambdaCDM(name=r'Flat $\Lambda$CDM', H0=h0*100. * u.km / u.s / u.Mpc, Om0=omega_matter)
        
        
        
        # Dictionary to save all observables in
        save_data = {}
        
        
        
#         stellar_masses = np.array([gal.masses['stellar'] for gal in obj.galaxies])
#         log_stellar_masses = np.log10(stellar_masses)
#         log_ssfr = np.log10(np.array([gal.sfr.to('Msun/Gyr')/gal.masses['stellar'].to('Msun') for gal in obj.galaxies]))
        
        stellar_masses = np.array([gal.masses['stellar'].to(mass_units) for gal in obj.galaxies if gal.central])
        log_stellar_masses = np.log10(stellar_masses)
        log_ssfr = np.log10(np.array([(gal.sfr/gal.masses['stellar']).to(ssfr_units) for gal in obj.galaxies if gal.central]))
        
        
        # Add eddington bias to calculated/measured/observed stellar masses (for GSMF only?)
        # See Kugel+23 (FLAMINGO: Calibrating large cosmological hydrodynamical simulations with machine learning)
        # Section 3.1.1 eqn 11
#         rng = np.random.default_rng()
#         mu_eddington = 0
#         sigma_eddington = min(0.070+0.071*redshift, 0.3) # dex
#         # print(sigma_eddington)
#         # eddington_bias = rng.lognormal(mean=mu_eddington, sigma=sigma_eddington)#, size=1)
#         eddington_bias = rng.normal(loc=mu_eddington, scale=sigma_eddington, size=len(log_stellar_masses))

#         # print([rng.normal(loc=mu_eddington, scale=sigma_eddington) for i in range(5)])

#         log_stellar_masses_eddington_biased = log_stellar_masses + eddington_bias
#         stellar_masses_eddington_biased = 10**log_stellar_masses_eddington_biased



        # sSFR Histograms
        x_units = ssfr_units
        y_units = unyt.dimensionless

        x_label = r'$\log(sSFR\,[Gyr^{-1}])$'
        y_label = r'fraction'
    
    
        logM_bins = np.array([9, 10, 11, 20])
        # log_ssfr_bins_ = np.arange(-3, 1, 0.2)
        log_ssfr_bins = np.arange(-3.1, 1, 0.2)
        log_ssfr_lim = -2.5 + 0.3*redshift

        log_ssfr_hists = gen.log_ssfr_hist_func(logM_bins, log_ssfr_bins, log_stellar_masses, log_ssfr, log_ssfr_lim=log_ssfr_lim)
        
        for ii in range(len(log_ssfr_hists)):
            save_data['ssfr_hist_bin_%s' % ii] = {
                'x':log_ssfr_hists[ii][0] * x_units,
                'xerr':np.zeros(len(log_ssfr_hists[ii][0])) * x_units,
                'y':log_ssfr_hists[ii][1] * y_units,
                'yerr':log_ssfr_hists[ii][2] * y_units,
                'x_label':x_label,
                'y_label':y_label,
                'name':'sSFR Histogram',
            }
            
            
        
        
        # Quenched fraction plotted as a function of stellar mass
        x_units = mass_units
        y_units = unyt.dimensionless

        x_label = r'$\log(M_{\star}\,/\,%s)$' % x_units.units.latex_representation()
        y_label = r'quenched fraction'
        
        
        logM_bins = np.array([9, 10, 11, 20])
        log_ssfr_bins = np.array([-np.inf, log_ssfr_lim])
    
        quenched_fraction = gen.ssfr_fraction_func(log_stellar_masses, log_ssfr, log_ssfr_bins, dlogM=1, 
                                                    min_logM=9, max_logM=13, calc_min_logM=False, calc_max_logM=False)
        quenched_fraction_v2 = gen.ssfr_fraction_func_v2(log_stellar_masses, log_ssfr, logM_bins, log_ssfr_bins)
        
        for ii in range(len(quenched_fraction)):
            save_data['quenched_fraction'] = {
                'x':quenched_fraction[ii][0] * x_units,
                'xerr':np.zeros(len(quenched_fraction[ii][0])) * x_units,
                'y':quenched_fraction[ii][1] * y_units,
                'yerr':quenched_fraction[ii][2] * y_units,
                'x_label':x_label,
                'y_label':y_label,
                'name':'Quenched Fraction',
            }
            
            save_data['quenched_fraction_v2'] = {
                'x':quenched_fraction_v2[ii][0] * x_units,
                'xerr':np.zeros(len(quenched_fraction_v2[ii][0])) * x_units,
                'y':quenched_fraction_v2[ii][1] * y_units,
                'yerr':quenched_fraction_v2[ii][2] * y_units,
                'x_label':x_label,
                'y_label':y_label,
                'name':'Quenched Fraction',
            }

        
        
        
        
        # Save to Velociraptor HDF5 files
        output_directory = SNAPDIR

        if not os.path.exists(output_directory):
            os.mkdir(output_directory)

        for key, val in save_data.items():
            output_filename = "%s_%s_%04d.hdf5" % (SIM, key, j)

            comment = f"h-corrected for SWIFT using cosmology: {cosmology.name}. "
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