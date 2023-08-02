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

from velociraptor.observations.objects import ObservationalData, MultiRedshiftObservationalData
#from astropy.cosmology import WMAP7 as cosmology
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

##########################################################################################################################
#### This function generates galaxy stellar mass functions (GSMF) from ceasar files for the given simulation in a given directory.
#### Call it using python gen_bhmsm.py SNAPDIR MODEL1 MODEL2...   XX
#### Call it using python gen_bhmsm.py MODELDIR SIMDIR SIM
#### Based on Renier Hough's "Create_caesarfiles.py"
##########################################################################################################################

parser = argparse.ArgumentParser()
#parser.add_argument('myArgs', nargs='*')
#parser.add_argument('myArgs', nargs=3)
parser.add_argument('--modeldir')
parser.add_argument('--simdir')
parser.add_argument('--sim', type=str)
args = parser.parse_args()
#MODELDIR = args.myArgs[0]
#NBINS = int(args.myArgs[1])
#SIM = []
#for i in range(1,len(args.myArgs)):
#    SIM.append(args.myArgs[i])
#SIMDIR = args.myArgs[1]
#SIM = args.myArgs[2]

MODELDIR = args.modeldir
SIMDIR = args.simdir
SIM = args.sim


#SNAPLIST = list(range(0,16))
SNAPLIST = [15]

print (SNAPLIST)

def gen_observable(SNAPLIST, SNAPDIR, SIM):
    print (SNAPDIR)
    '''
    SNAPLIST: list of the snapshot numbers
    SIM: Sim name, e.g. m50n512, used in snapshot name (assumes Simba naming convention)
    SNAPDIR = Directory where snapshots are located
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
        OUTFILE = os.path.join(SNAPDIR, '%s_bhmsm_%04d.txt' % (SIM,j))

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


        
        # Dictionary to save all observables in
        save_data = {}

        

        ## Simulation data for observable
        gal_stellar_masses = np.array([gal.masses['stellar'] for gal in obj.galaxies])    
        gal_bh_masses = np.array([gal.masses['bh'] for gal in obj.galaxies])

        min_log_gal_stellar_mass = 10
        max_log_gal_stellar_mass = 12

        #    gal_stellar_masses_for_observable = gal_stellar_masses[gal_bh_masses > 0]
        #    gal_bh_masses_for_observable = gal_bh_masses[gal_bh_masses > 0]

        #    gal_stellar_masses_for_observable = gal_stellar_masses_for_observable[gal_stellar_masses_for_observable >= 10**min_log_gal_stellar_mass]
        #    gal_bh_masses_for_observable = gal_bh_masses_for_observable[gal_stellar_masses_for_observable >= 10**min_log_gal_stellar_mass]

        gal_stellar_masses_no_zeros = gal_stellar_masses[gal_bh_masses > 0]
        gal_bh_masses_no_zeros = gal_bh_masses[gal_bh_masses > 0]

        gal_stellar_masses_for_observable = gal_stellar_masses_no_zeros[
            np.where(np.logical_and(
                gal_stellar_masses_no_zeros >= 10**min_log_gal_stellar_mass, 
                gal_stellar_masses_no_zeros <= 10**max_log_gal_stellar_mass))]
        gal_bh_masses_for_observable = gal_bh_masses_no_zeros[
            np.where(np.logical_and(
                gal_stellar_masses_no_zeros >= 10**min_log_gal_stellar_mass, 
                gal_stellar_masses_no_zeros <= 10**max_log_gal_stellar_mass))]
        
        log_gal_stellar_masses_for_observable = np.log10(gal_stellar_masses_for_observable)
        log_gal_bh_masses_for_observable = np.log10(gal_bh_masses_for_observable)



        save_data['data'] = {'x':gal_stellar_masses_for_observable * unyt.Solar_Mass,
                             'xerr':None,
                             'y':gal_bh_masses_for_observable * unyt.Solar_Mass,
                             'yerr':None}




        # Fit from KH13 (pg. 57, eqn 10)
        # They use a symmetric, least-squares ﬁt

        #    def Mbh_KH13_fit_fn(Mbulge, Mbh_Mbulge_ratio, power):
        #        Mbh_val = 1e9 * Mbh_Mbulge_ratio * (Mbulge/1e11)**power  # Msun
        #        return Mbh_val

        def log_Mbh_KH13_fit_fn(log_Mbulge, log_Mbh_Mbulge_ratio, power):
            log_Mbh_val = 9 + log_Mbh_Mbulge_ratio + power*(log_Mbulge - 11)  # log(M/Msun)
            return log_Mbh_val


        #    def Mbh_KH13_fit_fn_with_errors(Mbulge, params, param_errs):
        #        Mbh_val = 1e9 * params[0] * (Mbulge/1e11)**params[1]  # Msun
        #        Mbh_err = Mbh_val * np.sqrt( (param_errs[0]/params[0])**2 + ( np.log(Mbulge/1e11) * param_errs[1] )**2 ) # Msun
        #        return Mbh_val, Mbh_err

        def log_Mbh_KH13_fit_fn_with_errors(log_Mbulge, params, param_errs):
            log_Mbh_val = 9 + params[0] + params[1]*(log_Mbulge - 11)  # log(M/Msun)
            log_Mbh_err = np.sqrt( (param_errs[0])**2 + (param_errs[1]*(log_Mbulge - 11))**2 )
            return log_Mbh_val, log_Mbh_err


        #    popt, pcov = curve_fit(Mbh_KH13_fit_fn, gal_stellar_masses_for_observable, gal_bh_masses_for_observable)
        #    perr = np.sqrt(np.abs(np.diag(pcov)))
        popt, pcov = curve_fit(log_Mbh_KH13_fit_fn, np.log10(gal_stellar_masses_for_observable), np.log10(gal_bh_masses_for_observable))
        perr = np.sqrt(np.abs(np.diag(pcov)))

        N_points = 10
        gal_stellar_masses_fit = np.logspace(10, 12, N_points)
        log_gal_stellar_masses_fit = np.log10(gal_stellar_masses_fit)
        #    gal_bh_masses_fit, err_gal_bh_masses_fit = Mbh_KH13_fit_fn_with_errors(gal_stellar_masses_fit, popt, perr)
        log_gal_bh_masses_fit, log_err_gal_bh_masses_fit = log_Mbh_KH13_fit_fn_with_errors(log_gal_stellar_masses_fit, popt, perr)

        gal_bh_masses_fit = 10**log_gal_bh_masses_fit
        err_gal_bh_masses_fit = log_err_gal_bh_masses_fit * np.log(10) * log_gal_bh_masses_fit

        save_data['fit_sparse'] = {'x':gal_stellar_masses_fit * unyt.Solar_Mass,
                            'xerr':None,
                            'y':gal_bh_masses_fit * unyt.Solar_Mass,
                            'yerr':err_gal_bh_masses_fit * unyt.Solar_Mass}
        save_data['log_fit_sparse'] = {'x':log_gal_stellar_masses_fit * unyt.Solar_Mass,
                            'xerr':None,
                            'y':log_gal_bh_masses_fit * unyt.Solar_Mass,
                            'yerr':(10**log_gal_bh_masses_fit - 10**(log_gal_bh_masses_fit - log_err_gal_bh_masses_fit), 
                                    -10**log_gal_bh_masses_fit + 10**(log_gal_bh_masses_fit + log_err_gal_bh_masses_fit)) *
                                       unyt.Solar_Mass}



        
        
        
        # Bin data and find means/medians in each bin

        N_bins = 10
        #gal_stellar_masses_bins = np.logspace(min(log_gal_stellar_masses)-1e-5, max(log_gal_stellar_masses)+1e-5, N_bins+1)
        #gal_stellar_masses_bins = np.logspace(10, max(log_gal_stellar_masses)+1e-5, N_bins+1)
        gal_stellar_masses_bins = np.logspace(min_log_gal_stellar_mass, max_log_gal_stellar_mass, N_bins+1)
        #    print(gal_stellar_masses_bins)
        dlog_gal_stellar_masses = np.log10(gal_stellar_masses_bins[1:]) - np.log10(gal_stellar_masses_bins[:-1])
        log_gal_stellar_masses_bins_ax = np.log10(gal_stellar_masses_bins[:-1]) + dlog_gal_stellar_masses/2.
        gal_stellar_masses_bins_ax = 10**log_gal_stellar_masses_bins_ax


        gal_bh_masses_binned = [gal_bh_masses_for_observable[np.where((gal_stellar_masses_for_observable > low) & (gal_stellar_masses_for_observable <= high))] for low, high in zip(gal_stellar_masses_bins[:-1], gal_stellar_masses_bins[1:])]
        log_gal_bh_masses_binned = [np.log10(gal_bh_masses_vals) for gal_bh_masses_vals in gal_bh_masses_binned]

        gal_bh_masses_bin_means = np.array([np.nanmean(gal_bh_masses_vals) for gal_bh_masses_vals in gal_bh_masses_binned])
        gal_bh_masses_bin_stds = np.array([np.nanstd(gal_bh_masses_vals) for gal_bh_masses_vals in gal_bh_masses_binned])

        gal_bh_masses_bin_medians = np.array([np.nanmedian(gal_bh_masses_vals) for gal_bh_masses_vals in gal_bh_masses_binned])
        gal_bh_masses_bin_quantiles_lo = np.array([np.nanquantile(gal_bh_masses_vals, 0.16) for gal_bh_masses_vals in gal_bh_masses_binned])
        gal_bh_masses_bin_quantiles_hi = np.array([np.nanquantile(gal_bh_masses_vals, 0.84) for gal_bh_masses_vals in gal_bh_masses_binned])
        
        log_gal_bh_masses_bin_means = np.array([np.nanmean(log_gal_bh_masses_vals) for log_gal_bh_masses_vals in log_gal_bh_masses_binned])
        log_gal_bh_masses_bin_stds = np.array([np.nanstd(log_gal_bh_masses_vals) for log_gal_bh_masses_vals in log_gal_bh_masses_binned])

        log_gal_bh_masses_bin_medians = np.array([np.nanmedian(log_gal_bh_masses_vals) for log_gal_bh_masses_vals in log_gal_bh_masses_binned])
        log_gal_bh_masses_bin_quantiles_lo = np.array([np.nanquantile(log_gal_bh_masses_vals, 0.16) for log_gal_bh_masses_vals in log_gal_bh_masses_binned])
        log_gal_bh_masses_bin_quantiles_hi = np.array([np.nanquantile(log_gal_bh_masses_vals, 0.84) for log_gal_bh_masses_vals in log_gal_bh_masses_binned])


        save_data['binned_mean'] = {'x':gal_stellar_masses_bins_ax * unyt.Solar_Mass,
                                    'xerr':None,
                                    'y':gal_bh_masses_bin_means * unyt.Solar_Mass,
                                    'yerr':gal_bh_masses_bin_stds * unyt.Solar_Mass}
        save_data['binned_median'] = {'x':gal_stellar_masses_bins_ax * unyt.Solar_Mass,
                                      'xerr':None,
                                      'y':gal_bh_masses_bin_medians * unyt.Solar_Mass,
                                      'yerr':(gal_bh_masses_bin_medians-gal_bh_masses_bin_quantiles_lo,
                                              gal_bh_masses_bin_medians+gal_bh_masses_bin_quantiles_hi) * unyt.Solar_Mass}
        
        save_data['binned_log_mean'] = {'x':10**log_gal_stellar_masses_bins_ax * unyt.Solar_Mass,
                                    'xerr':None,
                                    'y':10**log_gal_bh_masses_bin_means * unyt.Solar_Mass,
                                    'yerr':(10**log_gal_bh_masses_bin_means - 
                                            10**(log_gal_bh_masses_bin_means - log_gal_bh_masses_bin_stds), 
                                            -10**log_gal_bh_masses_bin_means + 
                                            10**(log_gal_bh_masses_bin_means + log_gal_bh_masses_bin_stds)) * unyt.Solar_Mass}
        save_data['binned_log_median'] = {'x':10**log_gal_stellar_masses_bins_ax * unyt.Solar_Mass,
                                      'xerr':None,
                                      'y':10**log_gal_bh_masses_bin_medians * unyt.Solar_Mass,
                                      'yerr':(10**log_gal_bh_masses_bin_medians - 10**log_gal_bh_masses_bin_quantiles_lo, 
                                              -10**log_gal_bh_masses_bin_medians + 10**log_gal_bh_masses_bin_quantiles_hi) * 
                                          unyt.Solar_Mass}
        




        
        

        # Bin data and find mean/median in each bin, with overlapping bins (essentially rolling/running mean with fixed bins and window size,
        # instead of just fixed window in number of points)
        # Add in arbitrary window size in terms of x axis units (or fractions of bins?)

        N_bins = 150
        half_window_size_in_bins = 20
        #gal_stellar_masses_bins = np.logspace(9, 12+np.log10(2), N_bins+1)
        #gal_stellar_masses_bins = np.logspace(min(log_gal_stellar_masses)-1e-5, max(log_gal_stellar_masses)+1e-5, N_bins+1)
        #gal_stellar_masses_bins = np.logspace(10, max(log_gal_stellar_masses)+1e-5, N_bins+1)
        gal_stellar_masses_bins = np.logspace(min_log_gal_stellar_mass, min_log_gal_stellar_mass, N_bins+1)
        dlog_gal_stellar_masses = np.log10(gal_stellar_masses_bins[1:]) - np.log10(gal_stellar_masses_bins[:-1])
        log_gal_stellar_masses_bins_ax_windowed = np.log10(gal_stellar_masses_bins[half_window_size_in_bins:-half_window_size_in_bins-1]) + dlog_gal_stellar_masses[half_window_size_in_bins:-half_window_size_in_bins]/2.
        gal_stellar_masses_bins_ax_windowed = 10**log_gal_stellar_masses_bins_ax_windowed

        gal_bh_masses_binned = [gal_bh_masses_for_observable[np.where((gal_stellar_masses_for_observable > low) & (gal_stellar_masses_for_observable <= high))] for low, high in zip(gal_stellar_masses_bins[:-1], gal_stellar_masses_bins[1:])]

        gal_bh_masses_windowed = []
        for ii in range(half_window_size_in_bins, len(gal_bh_masses_binned)-half_window_size_in_bins):
            gal_bh_masses_windowed_curr = np.array([])
            for jj in range(-half_window_size_in_bins, half_window_size_in_bins+1): # check this range
                gal_bh_masses_windowed_curr = np.append(gal_bh_masses_windowed_curr, gal_bh_masses_binned[ii+jj])

            gal_bh_masses_windowed.append(gal_bh_masses_windowed_curr)

            
        gal_bh_masses_bin_means_windowed = np.array([np.nanmean(gal_bh_masses_vals) for gal_bh_masses_vals in gal_bh_masses_windowed])
        gal_bh_masses_bin_stds_windowed = np.array([np.nanstd(gal_bh_masses_vals) for gal_bh_masses_vals in gal_bh_masses_windowed])

        gal_bh_masses_bin_medians_windowed = np.array([np.nanmedian(gal_bh_masses_vals) for gal_bh_masses_vals in gal_bh_masses_windowed])
        gal_bh_masses_bin_quantiles_lo_windowed = np.array([np.nanquantile(gal_bh_masses_vals, 0.16) for gal_bh_masses_vals in gal_bh_masses_windowed])
        gal_bh_masses_bin_quantiles_hi_windowed = np.array([np.nanquantile(gal_bh_masses_vals, 0.84) for gal_bh_masses_vals in gal_bh_masses_windowed])


        save_data['windowed_mean'] = {'x':gal_stellar_masses_bins_ax_windowed * unyt.Solar_Mass,
                                      'xerr':None,
                                      'y':gal_bh_masses_bin_means_windowed * unyt.Solar_Mass,
                                      'yerr':gal_bh_masses_bin_stds_windowed * unyt.Solar_Mass}
        save_data['windowed_median'] = {'x':gal_stellar_masses_bins_ax_windowed * unyt.Solar_Mass,
                                        'xerr':None,
                                        'y':gal_bh_masses_bin_medians_windowed * unyt.Solar_Mass,
                                        'yerr':(gal_bh_masses_bin_medians_windowed-gal_bh_masses_bin_quantiles_lo_windowed,
                                                gal_bh_masses_bin_medians_windowed+gal_bh_masses_bin_quantiles_hi_windowed) * unyt.Solar_Mass}







        # Variable binning such that there are exact_countLY X objects (cannot be set) in N bins (must be set)
        #exact_count_N_obj = 3
        exact_count_N_bins = 5

        gal_stellar_masses_exact_count_binned, gal_stellar_masses_exact_count_bins = pd.qcut(gal_stellar_masses_for_observable, q=exact_count_N_bins, labels=None, retbins=True)

        dlog_gal_stellar_masses_exact_count_bins = np.log10(gal_stellar_masses_exact_count_bins[1:]) - np.log10(gal_stellar_masses_exact_count_bins[:-1])
        log_gal_stellar_masses_exact_count_bins_ax = np.log10(gal_stellar_masses_exact_count_bins[:-1]) + dlog_gal_stellar_masses_exact_count_bins/2.
        gal_stellar_masses_exact_count_bins_ax = 10**log_gal_stellar_masses_exact_count_bins_ax

        gal_bh_masses_exact_count_binned = [gal_bh_masses_for_observable[
            np.where((gal_stellar_masses_for_observable > low) & (gal_stellar_masses_for_observable <= high))]
                                      for low, high in zip(gal_stellar_masses_exact_count_bins[:-1], gal_stellar_masses_exact_count_bins[1:])]
        log_gal_bh_masses_exact_count_binned = [np.log10(gal_bh_masses_vals) for gal_bh_masses_vals in gal_bh_masses_exact_count_binned]

        
        gal_bh_masses_exact_count_bin_means = np.array([np.nanmean(gal_bh_masses_vals) for gal_bh_masses_vals in gal_bh_masses_exact_count_binned])
        gal_bh_masses_exact_count_bin_stds = np.array([np.nanstd(gal_bh_masses_vals) for gal_bh_masses_vals in gal_bh_masses_exact_count_binned])

        gal_bh_masses_exact_count_bin_medians = np.array([np.nanmedian(gal_bh_masses_vals) for gal_bh_masses_vals in gal_bh_masses_exact_count_binned])
        gal_bh_masses_exact_count_bin_quantiles_lo = np.array([np.nanquantile(gal_bh_masses_vals, 0.16) for gal_bh_masses_vals in gal_bh_masses_exact_count_binned])
        gal_bh_masses_exact_count_bin_quantiles_hi = np.array([np.nanquantile(gal_bh_masses_vals, 0.84) for gal_bh_masses_vals in gal_bh_masses_exact_count_binned])
        
        
        log_gal_bh_masses_exact_count_bin_means = np.array([np.nanmean(log_gal_bh_masses_vals) for log_gal_bh_masses_vals in log_gal_bh_masses_exact_count_binned])
        log_gal_bh_masses_exact_count_bin_stds = np.array([np.nanstd(log_gal_bh_masses_vals) for log_gal_bh_masses_vals in log_gal_bh_masses_exact_count_binned])

        log_gal_bh_masses_exact_count_bin_medians = np.array([np.nanmedian(log_gal_bh_masses_vals) for log_gal_bh_masses_vals in log_gal_bh_masses_exact_count_binned])
        log_gal_bh_masses_exact_count_bin_quantiles_lo = np.array([np.nanquantile(log_gal_bh_masses_vals, 0.16) for log_gal_bh_masses_vals in log_gal_bh_masses_exact_count_binned])
        log_gal_bh_masses_exact_count_bin_quantiles_hi = np.array([np.nanquantile(log_gal_bh_masses_vals, 0.84) for log_gal_bh_masses_vals in log_gal_bh_masses_exact_count_binned])
        

        save_data['exact_count_binned_mean'] = {'x':gal_stellar_masses_exact_count_bins_ax * unyt.Solar_Mass,
                                    'xerr':None,
                                    'y':gal_bh_masses_exact_count_bin_means * unyt.Solar_Mass,
                                    'yerr':gal_bh_masses_exact_count_bin_stds * unyt.Solar_Mass}
        save_data['exact_count_binned_median'] = {'x':gal_stellar_masses_exact_count_bins_ax * unyt.Solar_Mass,
                                      'xerr':None,
                                      'y':gal_bh_masses_exact_count_bin_medians * unyt.Solar_Mass,
                                      'yerr':(gal_bh_masses_exact_count_bin_medians-gal_bh_masses_exact_count_bin_quantiles_lo,
                                              gal_bh_masses_exact_count_bin_medians+gal_bh_masses_exact_count_bin_quantiles_hi) * 
                                                  unyt.Solar_Mass}
        
        save_data['exact_count_binned_log_mean'] = {'x':10**log_gal_stellar_masses_exact_count_bins_ax * unyt.Solar_Mass,
                                    'xerr':None,
                                    'y':10**log_gal_bh_masses_exact_count_bin_means * unyt.Solar_Mass,
                                    'yerr':(10**log_gal_bh_masses_exact_count_bin_means - 
                                            10**(log_gal_bh_masses_exact_count_bin_means - log_gal_bh_masses_exact_count_bin_stds), 
                                            -10**log_gal_bh_masses_exact_count_bin_means + 
                                            10**(log_gal_bh_masses_exact_count_bin_means + log_gal_bh_masses_exact_count_bin_stds)) *
                                                    unyt.Solar_Mass}
        save_data['exact_count_binned_log_median'] = {'x':10**log_gal_stellar_masses_exact_count_bins_ax * unyt.Solar_Mass,
                                      'xerr':None,
                                      'y':10**log_gal_bh_masses_exact_count_bin_medians * unyt.Solar_Mass,
                                      'yerr':(10**log_gal_bh_masses_exact_count_bin_medians - 
                                              10**log_gal_bh_masses_exact_count_bin_quantiles_lo,
                                              -10**log_gal_bh_masses_exact_count_bin_medians + 
                                              10**log_gal_bh_masses_exact_count_bin_quantiles_hi) * 
                                                  unyt.Solar_Mass}








        # Variable binning such that there are AT LEAST N objects in each bin
        min_N_obj = 3
        min_log_bin_width = 0.2
        log_delta_x = 0.02

        x_lo = 10**min_log_gal_stellar_mass
        #    x_lo = max(10**min_log_gal_stellar_mass, min(gal_stellar_masses_for_observable))
        x_hi = 10**max_log_gal_stellar_mass
        #    x_hi = min(10**max_log_gal_stellar_mass, max(gal_stellar_masses_for_observable))

        x = x_lo
        bin_edges = [x]

        hit_limit = False

        while True:
            x_up = deepcopy(x)
            while True:
                if (x_up == x):
                    log_x_up = np.log10(x_up) + min_log_bin_width + log_delta_x
                else:
                    log_x_up = np.log10(x_up) + log_delta_x
                x_up = 10**log_x_up
                #print("x:", x, " x_up:", x_up)

                log_bin_width = log_x_up - np.log10(x)
                #print("log_bin_width:", log_bin_width)

                gal_bh_masses_binned = gal_bh_masses_for_observable[np.where((gal_stellar_masses_for_observable > x) 
                                                                          & (gal_stellar_masses_for_observable <= x_up))]

                N_obj = len(gal_bh_masses_binned)
                #print("N_obj:", N_obj)


                if not hit_limit:
                    if (N_obj >= min_N_obj and log_bin_width >= min_log_bin_width):
                        #print("Conditions achieved")
                        break
                else:
                    #print("x_up one past x_hi")
                    break

                if (x_up > x_hi):
                    #print("x_up hit x_hi")
                    hit_limit = True

                #print("Conditions not achieved")
                #print()

            x = x_up
            if (x > x_hi and not hit_limit):
                #print("x hit x_hi")
                break
            elif (x > x_hi and hit_limit):
                #print("x hit x_hi")
                #print("Adding last upper bin edge")
                bin_edges.append(x)
                break

            #print("Adding upper bin edge")
            bin_edges.append(x)
            #print()
            #print()

        #print()
        #print()
        #print("Done")
        #print(bin_edges)


        gal_stellar_masses_min_count_bins = bin_edges

        dlog_gal_stellar_masses_min_count_bins = np.log10(gal_stellar_masses_min_count_bins[1:]) - np.log10(gal_stellar_masses_min_count_bins[:-1])
        log_gal_stellar_masses_min_count_bins_ax = np.log10(gal_stellar_masses_min_count_bins[:-1]) + dlog_gal_stellar_masses_min_count_bins/2.
        gal_stellar_masses_min_count_bins_ax = 10**log_gal_stellar_masses_min_count_bins_ax

        gal_bh_masses_min_count_binned = [gal_bh_masses_for_observable[
            np.where((gal_stellar_masses_for_observable > low) & (gal_stellar_masses_for_observable <= high))]
                                      for low, high in zip(gal_stellar_masses_min_count_bins[:-1], gal_stellar_masses_min_count_bins[1:])]
        log_gal_bh_masses_min_count_binned = [np.log10(gal_bh_masses_vals) for gal_bh_masses_vals in gal_bh_masses_min_count_binned]

        #print([len(bin) for bin in gal_bh_masses_min_count_binned])



        gal_bh_masses_min_count_bin_means = np.array([np.nanmean(gal_bh_masses_vals) for gal_bh_masses_vals in gal_bh_masses_min_count_binned])
        gal_bh_masses_min_count_bin_stds = np.array([np.nanstd(gal_bh_masses_vals) for gal_bh_masses_vals in gal_bh_masses_min_count_binned])

        gal_bh_masses_min_count_bin_medians = np.array([np.nanmedian(gal_bh_masses_vals) for gal_bh_masses_vals in gal_bh_masses_min_count_binned])

        gal_bh_masses_min_count_bin_quantiles_lo = np.array([np.nanquantile(gal_bh_masses_vals, 0.16) for gal_bh_masses_vals in gal_bh_masses_min_count_binned])
        gal_bh_masses_min_count_bin_quantiles_hi = np.array([np.nanquantile(gal_bh_masses_vals, 0.84) for gal_bh_masses_vals in gal_bh_masses_min_count_binned])
                                              
                                              
        log_gal_bh_masses_min_count_bin_means = np.array([np.nanmean(log_gal_bh_masses_vals) for log_gal_bh_masses_vals in log_gal_bh_masses_min_count_binned])
        log_gal_bh_masses_min_count_bin_stds = np.array([np.nanstd(log_gal_bh_masses_vals) for log_gal_bh_masses_vals in log_gal_bh_masses_min_count_binned])

        log_gal_bh_masses_min_count_bin_medians = np.array([np.nanmedian(log_gal_bh_masses_vals) for log_gal_bh_masses_vals in log_gal_bh_masses_min_count_binned])
        log_gal_bh_masses_min_count_bin_quantiles_lo = np.array([np.nanquantile(log_gal_bh_masses_vals, 0.16) for log_gal_bh_masses_vals in log_gal_bh_masses_min_count_binned])
        log_gal_bh_masses_min_count_bin_quantiles_hi = np.array([np.nanquantile(log_gal_bh_masses_vals, 0.84) for log_gal_bh_masses_vals in log_gal_bh_masses_min_count_binned])



        save_data['min_count_binned_mean'] = {'x':gal_stellar_masses_min_count_bins_ax * unyt.Solar_Mass,
                                    'xerr':None,
                                    'y':gal_bh_masses_min_count_bin_means * unyt.Solar_Mass,
                                    'yerr':gal_bh_masses_min_count_bin_stds * unyt.Solar_Mass}
        save_data['min_count_binned_median'] = {'x':gal_stellar_masses_min_count_bins_ax * unyt.Solar_Mass,
                                      'xerr':None,
                                      'y':gal_bh_masses_min_count_bin_medians * unyt.Solar_Mass,
                                      'yerr':(gal_bh_masses_min_count_bin_medians-gal_bh_masses_min_count_bin_quantiles_lo,
                                              gal_bh_masses_min_count_bin_medians+gal_bh_masses_min_count_bin_quantiles_hi) * unyt.Solar_Mass}
                                              
                                              
        save_data['min_count_binned_log_mean'] = {'x':10**log_gal_stellar_masses_min_count_bins_ax * unyt.Solar_Mass,
                                    'xerr':None,
                                    'y':10**log_gal_bh_masses_min_count_bin_means * unyt.Solar_Mass,
                                    'yerr':(10**log_gal_bh_masses_min_count_bin_means - 
                                            10**(log_gal_bh_masses_min_count_bin_means - log_gal_bh_masses_min_count_bin_stds), 
                                            -10**log_gal_bh_masses_min_count_bin_means + 
                                            10**(log_gal_bh_masses_min_count_bin_means + log_gal_bh_masses_min_count_bin_stds)) * 
                                                  unyt.Solar_Mass}
        save_data['min_count_binned_log_median'] = {'x':10**log_gal_stellar_masses_min_count_bins_ax * unyt.Solar_Mass,
                                      'xerr':None,
                                      'y':10**log_gal_bh_masses_min_count_bin_medians * unyt.Solar_Mass,
                                      'yerr':(10**log_gal_bh_masses_min_count_bin_medians - 
                                              10**log_gal_bh_masses_min_count_bin_quantiles_lo,
                                              -10**log_gal_bh_masses_min_count_bin_medians + 
                                              10**log_gal_bh_masses_min_count_bin_quantiles_hi) * unyt.Solar_Mass}







        # Save to HDF5 files
        output_directory = SNAPDIR

        if not os.path.exists(output_directory):
            os.mkdir(output_directory)

        for key, val in save_data.items():
            output_filename = "%s_bhmsm_%s_%04d.hdf5" % (SIM, key, j)

            comment = f"h-corrected for SWIFT using cosmology: {cosmology.name}."
            citation = ""
            bibcode = ""
            name = "BHMSM (%s) relation from swimba" % key
            plot_as = "points"
            redshift_lower = redshift
            redshift_upper = redshift
            h = cosmology.h


            x = val['x']
            xerr = val['xerr']

            y = val['y']
            yerr = val['yerr']

            processed = ObservationalData()
            processed.associate_x(x, scatter=xerr, comoving=True, description="Galaxy Stellar Mass")
            processed.associate_y(y, scatter=yerr, comoving=True, description="Galaxy Black Hole Mass")
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



'''
for i in range(0, len(SIM)):
#    gen_observable(SNAPLIST, SNAPDIR, SIM=SIM[i])    
    gen_observable(SNAPLIST, os.path.join(MODELDIR, SIMDIR[i]), SIM)
    
    print()
    print("DONE")
'''
                                              
gen_observable(SNAPLIST, os.path.join(MODELDIR, SIMDIR), SIM)

#print()
print("DONE")
print()
print()