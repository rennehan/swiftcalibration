USE_VERSION = 'v0.2b'
import os
import sys
import yt
import swiftsimio
import caesar
import unyt
import argparse
import numpy as np
import pandas as pd
from copy import copy, deepcopy

import gen_sim_data as gen

from velociraptor.observations.objects import ObservationalData, MultiRedshiftObservationalData
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

##########################################################################################################################
#### This function generates Cosmic Star Formation History (CSFH = Star Formation Rate Density (SFRD) across cosmic time
#### from caesar files for the given simulation in a given directory.
#### Call it using python gen_csfh.py --modeldir=path/to/model --simdir=SIMDIR --sim=SIM
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
sfr_output_units = unyt.msun / (unyt.year * unyt.Mpc**3)

x_units = unyt.dimensionless
y_units = unyt.msun / (unyt.year * unyt.Mpc**3)

x_label = r'z'
y_label = r'$\Psi=SFRD=\rho_{SFR}=\dot{\rho_{*}}\, [M_\odot\, yr^{-1}\, cMpc^{-3}]$'

x_label_log = r'$\log(1+z)$'
y_label_log = r'$\log(\Psi=SFRD=\rho_{SFR}=\dot{\rho_{*}}\, [M_\odot\, yr^{-1}\, cMpc^{-3}])$'




def load_sfr_data(simulation, preloaded=None):
    filename = f"{simulation}/SFR.txt"
    data = np.genfromtxt(filename).T

    print(filename)
    if preloaded is None:
        preloaded = {}
        preloaded['snapshot'] = "%s/%s_%04d.hdf5" % (simulation, SIM, SNAPLIST[0])
        preloaded['initial_snapshot'] = swiftsimio.load(preloaded['snapshot'])
        preloaded['units'] = preloaded['initial_snapshot'].units
        preloaded['boxsize'] = preloaded['initial_snapshot'].metadata.boxsize
        preloaded['box_volume'] = preloaded['boxsize'][0] * preloaded['boxsize'][1] * preloaded['boxsize'][2]

        preloaded['sfr_units'] = preloaded['initial_snapshot'].gas.star_formation_rates.units

    # a, Redshift, SFR
    return data[2], data[3], (data[7] * preloaded['sfr_units'] / preloaded['box_volume']).to(y_units)



def gen_observable(SNAPLIST, SNAPDIR, SIM):
    print (SNAPDIR)
    '''
    SNAPLIST: list of the snapshot numbers
    SNAPDIR = Directory where snapshots are located
    SIM: Sim name, e.g. m50n512, used in snapshot name (assumes Simba naming convention)
    '''
    
    
    
    # path to the file
    SNAP  = os.path.join(SNAPDIR, '%s_0000.hdf5' % SIM)
    if not os.path.exists(SNAP):
        print(SNAP, "does not exist")
        return

    CAESARFILE = os.path.join(SNAPDIR, '%s_caesar_%04d.hdf5' % (SIM, SNAPLIST[0]))

    ds = yt.load(SNAP)
    Vcom = ds.domain_width.in_units("Mpccm").prod()  # Comoving box volume in Mpc^3

    print('Loading caesar object')
    obj = caesar.load(CAESARFILE)
    redshift_caesar = obj.simulation.redshift
    h0_caesar = obj.simulation.hubble_constant
    boxsize_caesar = obj.simulation.boxsize.in_units("Mpccm")

    print("Comoving Volume=%s" % Vcom)
    print('z=%s' % redshift_caesar)
    print('h0=%s' % h0_caesar)


    # Cosmological Parameters
    omega_matter = ds.cosmology.omega_matter
    cosmology = FlatLambdaCDM(name=r'Flat $\Lambda$CDM', H0=h0_caesar*100. * u.km / u.s / u.Mpc, Om0=omega_matter)
    
    
    
    
    ## Actually calculate CSFH
    simulations_ = {"%s" % SNAPDIR: SIM}

    sim_keys = [k for k in simulations_.keys()]

    preloaded = None
    simulation_data = {k: load_sfr_data(k, preloaded) for k in simulations_.keys()}
    
    
    
    
    
    ## Dictionary in which to save all observables
    save_data = {}
    
    
    
    for simulation_ in simulation_data.keys():
        scale_factor_, redshift_, sfr_ = simulation_data[simulation_]
        name_ = simulations_[simulation_]
#         scale_factor = scale_factor_
#         redshift_ = redshift_
#         sfr_ = sfr_
    
    
    
        save_data['data'] = {
            'x':redshift_ * x_units,
            'xerr':np.zeros(len(redshift_)) * x_units,
            'y':sfr_.value * y_units,
            'yerr':np.zeros(len(sfr_.value)) * y_units,
            'x_label':x_label,
            'y_label':y_label,
            'plot_as':'points',
        }
        
        x_bin_centres, y_bin_means, y_bin_stds, y_bin_medians, y_bin_quantiles_lo, y_bin_quantiles_hi = gen.regular_bin(
            redshift_, sfr_.value, 
            dx=0.01, min_x=0, max_x=100,
            calc_min_x=False, calc_max_x=False)
        
        save_data['binned_median'] = {
            'x':x_bin_centres * x_units,
            'xerr':np.zeros(len(x_bin_centres)) * x_units,
            'y':y_bin_medians * y_units,
            'yerr':(np.abs(y_bin_medians - y_bin_quantiles_lo), np.abs(y_bin_medians - y_bin_quantiles_hi)) * y_units,
            'x_label':x_label,
            'y_label':y_label,
            'plot_as':'line',
        }
        
        save_data['binned_mean'] = {
            'x':x_bin_centres * x_units,
            'xerr':np.zeros(len(x_bin_centres)) * x_units,
            'y':y_bin_means * y_units,
            'yerr':y_bin_stds * y_units,
            'x_label':x_label,
            'y_label':y_label,
            'plot_as':'line',
        }
        
        
        save_data['log_data'] = {
            'x':np.log10(1+redshift_) * x_units,
            'xerr':np.zeros(len(np.log10(1+redshift_))) * x_units,
            'y':np.log10(sfr_) * y_units,
            'yerr':np.zeros(len(np.log10(sfr_))) * y_units,
            'x_label':x_label_log,
            'y_label':y_label_log,
            'plot_as':'points',
        }
        
        
        
        x_bin_centres, y_bin_means, y_bin_stds, y_bin_medians, y_bin_quantiles_lo, y_bin_quantiles_hi = gen.regular_bin(
            np.log10(1+redshift_), np.log10(sfr_), 
            dx=0.01, min_x=0, max_x=100,
            calc_min_x=False, calc_max_x=False)
        
        save_data['binned_log_median'] = {
            'x':x_bin_centres * x_units,
            'xerr':np.zeros(len(x_bin_centres)) * x_units,
            'y':y_bin_medians * y_units,
            'yerr':(np.abs(y_bin_medians - y_bin_quantiles_lo), np.abs(y_bin_medians - y_bin_quantiles_hi)) * y_units,
            'x_label':x_label,
            'y_label':y_label,
            'plot_as':'line',
        }
        
        save_data['binned_log_mean'] = {
            'x':x_bin_centres * x_units,
            'xerr':np.zeros(len(x_bin_centres)) * x_units,
            'y':y_bin_means * y_units,
            'yerr':y_bin_stds * y_units,
            'x_label':x_label,
            'y_label':y_label,
            'plot_as':'line',
        }
        
        
        
        # Save to HDF5 files
        output_directory = SNAPDIR

        if not os.path.exists(output_directory):
            os.mkdir(output_directory)

        for key, val in save_data.items():
            output_filename = "%s_csfh_%s.hdf5" % (name_, key)

            comment = f"h-corrected for SWIFT using cosmology: {cosmology}."
            citation = ""
            bibcode = ""
            name = "CSFH"
            plot_as = val['plot_as']
            redshift_lower = 0
            redshift_upper = 0


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
            processed.associate_redshift(0, redshift_lower, redshift_upper)
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