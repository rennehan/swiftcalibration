# Functions for generating observables from simulation data
import numpy as np
import pandas as pd
from copy import copy, deepcopy
from scipy.optimize import curve_fit


## GSMF Functions

def calc_bins(x, dx=0.5, min_x=0, max_x=100, calc_min_x=True, calc_max_x=True):
    if calc_min_x:
        min_x = np.nanmin(x)
    if calc_max_x:
        max_x = np.nanmax(x)
    
    # Makes sure the upper limit max_val is also included and not excluded
    max_bin = min_x + dx + dx * np.ceil((max_x - min_x)/float(dx))
    bins = np.arange(min_x, max_bin, dx)
    nbins = len(bins) - 1
    
    bin_centres = bins[0:-1] + dx/2.

    return bins, bin_centres, nbins

def calc_exact_count_bins(x, Nbins=5, min_x=0, max_x=100, calc_min_x=True, calc_max_x=True):
    if calc_min_x:
        min_x = np.nanmin(x)
    if calc_max_x:
        max_x = np.nanmax(x)
    
    x_ = x[np.where(np.logical_and(x >= min_x, x <= max_x))]
    
    x_binned, x_bins = pd.qcut(x_, q=Nbins, labels=None, retbins=True)
    dx_bins = x_bins[1:] - x_bins[:-1]
    x_bin_centres = x_bins[:-1] + dx_bins/2.
    
    return x_, x_binned, x_bins, x_bin_centres

def calc_min_count_bins(x, x_lo, x_hi, calc_min_x=True, calc_max_x=True,  min_N_obj=3, min_bin_width=0.4, delta_x=0.02):
    if calc_min_x:
        x_lo = np.nanmin(x)
    if calc_max_x:
        x_hi = np.nanmax(x)
    
    x_curr = x_lo
    x_bins = [x_curr]

    hit_limit = False

    while True:
        x_up = deepcopy(x_curr)
        while True:
            if (x_up == x_curr):
                x_up += min_bin_width
            else:
                x_up += delta_x
            print("x_curr:", x_curr, " x_up:", x_up)

            bin_width = x_up - x_curr
            print("bin width:", bin_width)

            x_binned_curr = x[np.where((x > x_curr) & (x <= x_up))]

            N_obj = len(x_binned_curr)
            print("N_obj:", N_obj)


            if not hit_limit:
                if (N_obj >= min_N_obj and bin_width >= min_bin_width):
                    print("Conditions achieved")
                    break
            else:
                print("x_up one past x_hi")
                break

            if (x_up > x_hi):
                print("x_up hit x_hi")
                hit_limit = True
#                 break

            print("Conditions not achieved")
            print()

        x_curr = x_up
        if x_curr > x_hi:
            print("x_curr over x_hi")
            print("Adding last upper bin edge")
            x_bins.append(x_curr)
            break

        print("Adding upper bin edge")
        x_bins.append(x_curr)
        print()
        print()

    print()
    print()
    print("Done")
    print("x bins:")
    print(x_bins)
    
    x_bins = np.array(x_bins)    
    dx_bins = x_bins[1:] - x_bins[:-1]
    x_bin_centres = x_bins[:-1] + dx_bins/2.
    
    return x_bins, x_bin_centres

    

def mass_function(mass, vol_com_Mpc, dlogM=0.5, min_logM=6, max_logM=14, calc_min_logM=True, calc_max_logM=True):
    logbins, logbin_centres, nbins = calc_bins(np.log10(mass), dx=dlogM, 
                                               min_x=min_logM, max_x=max_logM, 
                                               calc_min_x=calc_min_logM, calc_max_x=calc_max_logM)
    bins = 10**logbins
    bin_centres = 10**logbin_centres
    
    x = bin_centres    
    hist, edg = np.histogram(mass, bins=bins, range=(bins.min(), bins.max()))
    y = hist / (vol_com_Mpc * dlogM)
    return x, y, bins, nbins


def cosmic_variance(mass, pos, boxsize, vol_com_Mpc, dlogM=0.5, min_logM=6, max_logM=14, calc_min_logM=True, calc_max_logM=True):
    pos = np.floor(pos / (0.5 * boxsize)).astype(np.int)
    gal_index = pos[:, 0] + pos[:, 1] * 2 + pos[:, 2] * 4
    x, y, bins, nbins = mass_function(mass, vol_com_Mpc, dlogM=dlogM, min_logM=min_logM, max_logM=max_logM, 
                                      calc_min_logM=calc_min_logM, calc_max_logM=calc_max_logM)

    store_mf = np.zeros((nbins, 8))
    for i0 in range(8):
        m_s_ = mass[gal_index==i0]
        if len(m_s_) < 1: 
            continue
            
        if calc_min_logM and calc_max_logM:
            _x, _y, _bins, _nbins = mass_function(m_s_, vol_com_Mpc, dlogM=dlogM, 
                                                  min_logM=np.log10(np.nanmin(mass)), max_logM=np.log10(np.nanmax(mass)),
                                                  calc_min_logM=False, calc_max_logM=False)
        elif calc_min_logM and not calc_max_logM:
            _x, _y, _bins, _nbins = mass_function(m_s_, vol_com_Mpc, dlogM=dlogM, 
                                                  min_logM=np.log10(np.nanmin(mass)), max_logM=max_logM,
                                                  calc_min_logM=False, calc_max_logM=False)
        elif not calc_min_logM and calc_max_logM:
            _x, _y, _bins, _nbins = mass_function(m_s_, vol_com_Mpc, dlogM=dlogM, 
                                                  min_logM=min_logM, max_logM=np.log10(np.nanmax(mass)),
                                                  calc_min_logM=False, calc_max_logM=False)
        else:
            _x, _y, _bins, _nbins = mass_function(m_s_, vol_com_Mpc, dlogM=dlogM, 
                                                  min_logM=min_logM, max_logM=max_logM,
                                                  calc_min_logM=False, calc_max_logM=False)

        phi_  = np.log10(8 * _y)
        phi_  = np.where(phi_ < -100, np.log10(y), phi_)
        store_mf[:,i0] = phi_
    store_mf = np.ma.masked_invalid(store_mf)
    var = np.ma.std(store_mf, axis=1)
    return np.log10(x), np.log10(y), np.array(var)


def poisson_error(mass, vol_com_Mpc, dlogM=0.5, min_logM=6, max_logM=14, calc_min_logM=True, calc_max_logM=True):
    x, Phi, bins, nbins = mass_function(mass, vol_com_Mpc, dlogM=dlogM, min_logM=min_logM, max_logM=max_logM,
                                        calc_min_logM=calc_min_logM, calc_max_logM=calc_max_logM)
    Phi_err = np.sqrt(Phi * (vol_com_Mpc * dlogM))/(vol_com_Mpc * dlogM)
    
    logPhi_lo_err = np.abs(np.log10(Phi) - np.log10(Phi - Phi_err))
    logPhi_hi_err = np.abs(np.log10(Phi) - np.log10(Phi + Phi_err))
    
    logPhi_err = Phi_err/(np.log(10)*Phi)
    
    logPhi_err_v2 = np.log10(Phi) * Phi_err/Phi
    
    return logPhi_lo_err, logPhi_hi_err, logPhi_err, logPhi_err_v2


def mass_function_with_error(mass, pos, boxsize, vol_com_Mpc, dlogM=0.5, min_logM=6, max_logM=14, calc_min_logM=True, calc_max_logM=True):
    logx, logPhi, log_Phi_cv_err = cosmic_variance(mass, pos, boxsize, vol_com_Mpc, dlogM=dlogM, min_logM=min_logM, max_logM=max_logM,
                                                   calc_min_logM=calc_min_logM, calc_max_logM=calc_max_logM)
    
    logPhi_lo_poisson_err, logPhi_hi_poisson_err, logPhi_poisson_err, logPhi_poisson_err_v2 = poisson_error(
        mass, vol_com_Mpc, dlogM=dlogM, 
        min_logM=min_logM, max_logM=max_logM,
        calc_min_logM=calc_min_logM, 
        calc_max_logM=calc_max_logM)
    
    logPhi_total_lo_err = np.sqrt(log_Phi_cv_err**2 + logPhi_lo_poisson_err**2)
    logPhi_total_hi_err = np.sqrt(log_Phi_cv_err**2 + logPhi_hi_poisson_err**2)
    logPhi_total_err = np.sqrt(log_Phi_cv_err**2 + logPhi_poisson_err**2)
    logPhi_total_err_v2 = np.sqrt(log_Phi_cv_err**2 + logPhi_poisson_err_v2**2)
    
    return logx, logPhi, logPhi_total_lo_err, logPhi_total_hi_err, logPhi_total_err, logPhi_total_err_v2


def baldry12_log_gsmf_schecter_fit(log_x):
    # Baldry+12 eqn 6
    # values/parameters changed so 
    log_Mstar = 10.66  # +/- 0.05
    phi1 = 3.96e-3  # +/- 0.34e-3
    phi2 = 0.79e-3  # +/- 0.23e-3
    alpha1 = -0.35  # +/- 0.18
    alpha2 = -1.47  # +/- 0.05
            
    Delta_logM_logMstar = log_x - log_Mstar
    
    PhiM = np.log(10) * np.exp(-10**Delta_logM_logMstar) * 10**Delta_logM_logMstar * (
        phi1 * 10**(Delta_logM_logMstar * alpha1) + phi2 * 10**(Delta_logM_logMstar * alpha2))

    return np.log10(PhiM)




## sSFR Functions

# Could add in bin calculator like in ssfr_fraction_func
def log_ssfr_hist_func(logM_bins, log_ssfr_bins, logM, log_ssfr, log_ssfr_lim=-2.5):
    # sets all "quenched" galaxies into quenched bin
#     log_ssfr = np.maximum(log_ssfr, log_ssfr_lim)
    log_ssfr = np.where(log_ssfr <= log_ssfr_lim, log_ssfr_lim - 1e-10*np.abs(log_ssfr_lim), log_ssfr)
    
    log_ssfr_hists = []
    for ii in range(len(logM_bins)-1):
        log_ssfr_curr = log_ssfr[np.logical_and(logM>=logM_bins[ii], logM<logM_bins[ii+1])]
        counts, bin_edges = np.histogram(log_ssfr_curr, bins=log_ssfr_bins)
#         counts, bin_edges = np.histogram(log_ssfr_curr, bins=log_ssfr_bins, density=True)
#         counts, bin_edges = np.histogram(log_ssfr_curr, bins=log_ssfr_bins, normed=True)
        bin_centres = 0.5*(bin_edges[:-1]+bin_edges[1:])
        poisson_err = np.sqrt(counts)
        
        frac = counts/len(log_ssfr_curr)
#         frac = counts/np.sum(counts)
        poisson_err /= len(log_ssfr_curr)
#         poisson_err /= np.sum(counts)


        log_ssfr_hists.append([bin_centres, frac, poisson_err])
    
    return log_ssfr_hists


def ssfr_fraction_func(logM, log_ssfr, log_ssfr_bins, dlogM=0.5, min_logM=9, max_logM=14, calc_min_logM=True, calc_max_logM=True):
    logM_bins, logM_bin_centres, nbins = calc_bins(logM, dx=dlogM, 
                                                       min_x=min_logM, max_x=max_logM, 
                                                       calc_min_x=calc_min_logM, calc_max_x=calc_max_logM)
    
    result = []
    for ii in range(len(log_ssfr_bins)-1):
        ssfr_fraction = []
        ssfr_fraction_err = []
        for low, high in zip(logM_bins[:-1], logM_bins[1:]):
            log_ssfr_vals = log_ssfr[np.logical_and(logM>=low, logM<high)]

            # Sanitize values of -infs
    #         log_ssfr_vals = log_ssfr_vals[np.isfinite(log_ssfr_vals)]

            try:
                log_ssfr_vals_bin = log_ssfr_vals[np.logical_and(log_ssfr_vals>log_ssfr_bins[ii], 
                                                                      log_ssfr_vals<=log_ssfr_bins[ii+1])]
                ssfr_fraction_val_bin = len(log_ssfr_vals_bin)/len(log_ssfr_vals)

#                 poisson_err_bin = np.sqrt(len(log_ssfr_vals_bin))/len(log_ssfr_vals)
                poisson_err_bin = ssfr_fraction_val_bin * np.sqrt(len(log_ssfr_vals_bin)**(-1) + len(log_ssfr_vals)**(-1))
            except:
                ssfr_fraction_val_bin = 0
                poisson_err_bin = 0

            ssfr_fraction.append(ssfr_fraction_val_bin)
            ssfr_fraction_err.append(poisson_err_bin)
        
        result.append([logM_bin_centres, np.array(ssfr_fraction), np.array(ssfr_fraction_err)])
        
    return result



def ssfr_fraction_func_v2(logM, log_ssfr, logM_bins, log_ssfr_bins):
    logM_bin_centres = 0.5*(logM_bins[1:] + logM_bins[:-1])
    
    result = []
    for ii in range(len(log_ssfr_bins)-1):
        ssfr_fraction = []
        ssfr_fraction_err = []
        for low, high in zip(logM_bins[:-1], logM_bins[1:]):
            log_ssfr_vals = log_ssfr[np.logical_and(logM>=low, logM<high)]

            # Sanitize values of -infs
    #         log_ssfr_vals = log_ssfr_vals[np.isfinite(log_ssfr_vals)]

            try:
                log_ssfr_vals_bin = log_ssfr_vals[np.logical_and(log_ssfr_vals>log_ssfr_bins[ii], 
                                                                 log_ssfr_vals<=log_ssfr_bins[ii+1])]
                ssfr_fraction_val_bin = len(log_ssfr_vals_bin)/len(log_ssfr_vals)

#                 poisson_err_bin = np.sqrt(len(log_ssfr_vals_bin))/len(log_ssfr_vals)
                poisson_err_bin = ssfr_fraction_val_bin * np.sqrt(len(log_ssfr_vals_bin)**(-1) + len(log_ssfr_vals)**(-1))
            except:
                ssfr_fraction_val_bin = 0
                poisson_err_bin = 0

            ssfr_fraction.append(ssfr_fraction_val_bin)
            ssfr_fraction_err.append(poisson_err_bin)
        
        result.append([logM_bin_centres, np.array(ssfr_fraction), np.array(ssfr_fraction_err)])
        
    return result




## Scaling relationship binning functions

def regular_bin(x, y, dx=0.5, min_x=0, max_x=100, calc_min_x=True, calc_max_x=True):
    # Bin data and find means/medians + std devs/quantiles in each bin
    
    x_bins, x_bin_centres, x_nbins = calc_bins(x, dx=dx, min_x=min_x, max_x=max_x, 
                                               calc_min_x=calc_min_x, calc_max_x=calc_max_x)
    
    y_binned = [y[np.where((x >= low) & (x < high))] for low, high in zip(x_bins[:-1], x_bins[1:])]
    
    y_bin_means = np.array([np.nanmean(y_vals) for y_vals in y_binned])
    y_bin_stds = np.array([np.nanstd(y_vals) for y_vals in y_binned])

    y_bin_medians = np.array([np.nanmedian(y_vals) for y_vals in y_binned])
    y_bin_quantiles_lo = np.array([np.nanquantile(y_vals, 0.16) for y_vals in y_binned])
    y_bin_quantiles_hi = np.array([np.nanquantile(y_vals, 0.84) for y_vals in y_binned])
    
    return x_bin_centres, y_bin_means, y_bin_stds, y_bin_medians, y_bin_quantiles_lo, y_bin_quantiles_hi

def exact_count_bin(x, y, Nbins=5, min_x=0, max_x=100, calc_min_x=True, calc_max_x=True):
    # Variable binning such that there are EXACTLY X objects (cannot be set) in N bins (must be set)
    
    x_, x_binned, x_bins, x_bin_centres = calc_exact_count_bins(x, Nbins=Nbins, min_x=min_x, max_x=max_x, 
                                                                calc_min_x=calc_min_x, calc_max_x=calc_max_x)
    
    y_binned = [y[np.where((x >= low) & (x < high))] for low, high in zip(x_bins[:-1], x_bins[1:])]
    
    y_bin_means = np.array([np.nanmean(y_vals) for y_vals in y_binned])
    y_bin_stds = np.array([np.nanstd(y_vals) for y_vals in y_binned])

    y_bin_medians = np.array([np.nanmedian(y_vals) for y_vals in y_binned])
    y_bin_quantiles_lo = np.array([np.nanquantile(y_vals, 0.16) for y_vals in y_binned])
    y_bin_quantiles_hi = np.array([np.nanquantile(y_vals, 0.84) for y_vals in y_binned])
    
    return x_bin_centres, y_bin_means, y_bin_stds, y_bin_medians, y_bin_quantiles_lo, y_bin_quantiles_hi

def min_count_bin(x, y, x_lo, x_hi, calc_min_x=True, calc_max_x=True, min_N_obj=3, min_bin_width=0.4, delta_x=0.02):
    # Variable binning such that there are AT LEAST N objects in each bin
    # with a minimum bin width required so bins are not too narrow if there are many objects clustered together
    
    x_bins, x_bin_centres = calc_min_count_bins(x, x_lo, x_hi, calc_min_x=calc_min_x, calc_max_x=calc_max_x,
                                                min_N_obj=min_N_obj, min_bin_width=min_bin_width, delta_x=delta_x)    

    y_binned = [y[np.where((x >= low) & (x < high))] for low, high in zip(x_bins[:-1], x_bins[1:])]
    
    y_bin_means = np.array([np.nanmean(y_vals) for y_vals in y_binned])
    y_bin_stds = np.array([np.nanstd(y_vals) for y_vals in y_binned])

    y_bin_medians = np.array([np.nanmedian(y_vals) for y_vals in y_binned])
    y_bin_quantiles_lo = np.array([np.nanquantile(y_vals, 0.16) for y_vals in y_binned])
    y_bin_quantiles_hi = np.array([np.nanquantile(y_vals, 0.84) for y_vals in y_binned])
    
    return x_bin_centres, y_bin_means, y_bin_stds, y_bin_medians, y_bin_quantiles_lo, y_bin_quantiles_hi


def kh13_log_fit_fn(log_x, ratio, power):
    log_y = 9 + ratio + power*(log_x - 11)  # log(M/Msun)
    return log_y

def kh13_log_fit_fn_with_errors(log_x, params, param_errs):
    log_y = 9 + params[0] + params[1]*(log_x - 11)  # log(M/Msun)
    
    if param_errs.ndim==1:
        log_y_err = np.sqrt( (param_errs[0])**2 + (param_errs[1]*(log_x - 11))**2 )
    elif param_errs.ndim==2:
        log_y_err_lo = np.sqrt( (param_errs[0][0])**2 + (param_errs[1][0]*(log_x - 11))**2 )
        log_y_err_hi = np.sqrt( (param_errs[0][1])**2 + (param_errs[1][1]*(log_x - 11))**2 )
        log_y_err = [log_y_err_lo, log_y_err_hi]
    else:
        print('Dimensionality of param_errs is too high')
        return
    
    return log_y, log_y_err

def b18_log_fit_fn(log_x, a, b, c):
    log_y = a + b*(log_x - c)  # log(M/Msun)
    return log_y

def b18_log_fit_fn_with_errors(log_x, params, param_errs):
    log_y = params[0] + params[1]*(log_x - params[2])  # log(M/Msun)
    
    if param_errs.ndim==1:
        log_y_err = np.sqrt( (param_errs[0])**2 + (param_errs[1]*(log_x - params[2]))**2 )
    elif param_errs.ndim==2:
        log_y_err_lo = np.sqrt( (param_errs[0][0])**2 + (param_errs[1][0]*(log_x - params[2][0]))**2 )
        log_y_err_hi = np.sqrt( (param_errs[0][1])**2 + (param_errs[1][1]*(log_x - params[2][1]))**2 )
        log_y_err = [log_y_err_lo, log_y_err_hi]
    else:
        print('Dimensionality of param_errs is too high')
        return
    
    return log_y, log_y_err

def g23_log_fit_fn(log_x, a, b, c):
    log_y = a + b*(log_x - c)  # log(M/Msun)   
    return log_y

def g23_log_fit_fn_with_errors(log_x, params, param_errs):
    log_y = params[0] + params[1]*(log_x - params[2])  # log(M/Msun)
    
    if param_errs.ndim==1:
        log_y_err = np.sqrt( (param_errs[0])**2 + (param_errs[1]*(log_x - params[2]))**2 )
    elif param_errs.ndim==2:
        log_y_err_lo = np.sqrt( (param_errs[0][0])**2 + (param_errs[1][0]*(log_x - params[2][0]))**2 )
        log_y_err_hi = np.sqrt( (param_errs[0][1])**2 + (param_errs[1][1]*(log_x - params[2][1]))**2 )
        log_y_err = [log_y_err_lo, log_y_err_hi]
    else:
        print('Dimensionality of param_errs is too high')
        return    
    
    return log_y, log_y_err

def log_fit_fn(log_x, log_y, func=None, paper='kh13', determine_func=True):
    # For actually fitting
    if determine_func:
        if paper=='kh13':
            func = kh13_log_fit_fn
        elif paper=='b18':
            func = b18_log_fit_fn
        elif paper=='g23':
            func = g23_log_fit_fn
        else:
            print('Invalid paper')
            return
    
    popt, pcov = curve_fit(func, log_x, log_y)
    perr = np.sqrt(np.abs(np.diag(pcov)))
    
    return popt, perr

def log_fit_fn_with_errors(log_x_arr, params, param_errs, func=None, paper='kh13', determine_func=True):
    if determine_func:
        if paper=='kh13':
            func = kh13_log_fit_fn_with_errors
        elif paper=='b18':
            func = b18_log_fit_fn_with_errors
        elif paper=='g23':
            func = g23_log_fit_fn_with_errors
        else:
            print('Invalid paper')
            return
        
    log_y, log_y_err = func(log_x_arr, params, param_errs)
    
    return log_y, log_y_err

def log_fit_fn_with_errors_for_fitting(log_x_arr, log_x, log_y, func=None, paper='kh13', determine_func=True):
    popt, perr = log_fit_fn(log_x, log_y, func=func, paper=paper, determine_func=determine_func)
    
    log_y, log_y_err = log_fit_fn_with_errors(log_x_arr, popt, perr, func=func, paper=paper, determine_func=determine_func)
    
    return log_x_arr, log_y, log_y_err