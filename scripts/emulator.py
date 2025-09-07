import numpy as np
import matplotlib.pyplot as plt
import unyt

import os
from glob import glob
from pathlib import Path

# For Swift Emulator training and validation
import swiftemulator as se
from swiftemulator.design import latin
from swiftemulator.io.swift import load_parameter_files
from swiftemulator.emulators import gaussian_process, gaussian_process_bins, gaussian_process_mcmc, gaussian_process_one_dim, linear_model, multi_gaussian_process
from swiftemulator.mean_models.polynomial import PolynomialMeanModel
from swiftemulator.sensitivity import basic, cross_check, cross_check_bins
from swiftemulator.mocking import mock_sweep


# For reading in simulation observables as Velociraptor HDF5 files
from velociraptor.observations import load_observation, load_observations


# For saving emulator for later use
import dill


#snap_prefix = 'simba_s25n256'  # eg. 'simba_s12.5n128'
snap_prefix = 'snapshot'

# Main directory where calibration runs have been performed - Adjust this for whatever I am using
#root_dir = '/scratch/b/babul/stlock/swimba_s25n256_ps2020_32x'
root_dir =  '/scratch/b/babul/stlock/snapshots/HyenasC/L0/halo_3224_finished_runs'
# entropy_venv_4stlock@nia-login03:/scratch/b/babul/stlock/swimba_s25n256_ps2020_32x/calibrations/cali_0000$ 

# Subdirectory that holds actual calibration runs
# cali_dir = 'calibrations'
cali_dir = 'BAL_variations'

# Name of simulation observable, eg. 'gsmf'
obs_name = '151_entropy_profile_common'


sim_paths = [Path(x) for x in glob(os.path.join(root_dir, cali_dir, 'run*', snap_prefix+'_'+obs_name+'.hdf5')) 
             if os.path.exists(x) and 'run06' not in x]

cali_dirs = [Path(x) for x in glob(os.path.join(root_dir, cali_dir, 'run*')) 
             if os.path.exists(os.path.join(x, snap_prefix+'_'+obs_name+'.hdf5')) and 'run06' not in x]

cali_dirnames = {dirname.stem: dirname for dirname in cali_dirs}
cali_list = [dirname.stem for dirname in cali_dirs]
print(sim_paths)
print(cali_list)
print(len(cali_list))
print(sim_paths)
print(cali_dirs)

sim_dict = {}
thin = 1

for cali, sim_path in zip(cali_list, sim_paths):
    print(cali)
    sim_file = sim_path
    sim_data = load_observations(sim_file)[0]

    try:
        sim_data = load_observations(sim_file)[0]
    except:
        print('Failed')
        continue

    # Load simulation observable data
    
    x_ = sim_data.x.value
    y_ = sim_data.y.value
    y_err_ = sim_data.y_scatter.value

    
    # Sanitize data of -inf, +inf, NaN values
    # Data for different simulations may have different array lengths
    isfinite = np.isfinite(x_)
    isfinite_err = np.resize(isfinite, np.shape(y_err_))
    x_ = x_[isfinite]
    y_ = y_[isfinite]
    y_err_copy = y_err_[isfinite_err]
    if (y_err_copy.ndim != y_err_.ndim):
        y_err_ = y_err_copy.reshape(y_err_.ndim, -1)
    else:
        y_err_ = y_err_copy
    

    isfinite = np.isfinite(y_)
    isfinite_err = np.resize(isfinite, np.shape(y_err_))
    x_ = x_[isfinite]
    y_ = y_[isfinite]
    y_err_copy = y_err_[isfinite_err]
    if (y_err_copy.ndim != y_err_.ndim):
        y_err_ = y_err_copy.reshape(y_err_.ndim, -1)
    else:
        y_err_ = y_err_copy
    
    
    # Get 1D errors
    y_err_arr = np.array(y_err_)
    if (y_err_arr.ndim!=1):
        mean_y_err_ = np.nanmean(y_err_arr, axis=0)
        max_y_err_ = np.nanmax(y_err_arr, axis=0)
    else:
        mean_y_err_ = y_err_arr
        max_y_err_ = y_err_arr
    
    
    sim_dict[cali] = {'x':x_[::thin], 'y':y_[::thin],
                      'y_err':y_err_, 'mean_y_err':mean_y_err_[::thin], 'max_y_err':max_y_err_[::thin], 
                      'sim_data':sim_data}


sim_data = load_observations(sim_file)[0]
sim_info = {
    'X_LABEL':sim_data.x_description,
    'Y_LABEL':sim_data.y_description,
    'name':sim_data.name,
    'x_units':sim_data.x_units,
    'y_units':sim_data.y_units,
    'comment':sim_data.comment,
    'z':sim_data.redshift,
    'z_lo':sim_data.redshift_lower,
    'z_hi':sim_data.redshift_upper,
    'plot_as':sim_data.plot_as,
}

# Plot the observable to check it worked

thin_ = 1
# thin_ = 10 # for CSFH

for cali in sim_dict:
    print("x: ", sim_dict[cali]['x'])
    print("y: ", sim_dict[cali]['y'])
    print("yerr: ", sim_dict[cali]['max_y_err'])    
    
    plt.clf()
    plt.xlabel(sim_info['X_LABEL'])
    plt.ylabel(sim_info['Y_LABEL'])
    plt.title(r'%s, %s' % (cali, snap_prefix))
    plt.grid()
    plt.minorticks_on()
    plt.errorbar(sim_dict[cali]['x'][::thin_], sim_dict[cali]['y'][::thin_], yerr=sim_dict[cali]['max_y_err'][::thin_], 
                 fmt='o-', markersize=2)
    # plt.scatter(sim_dict[cali]['x'][::thin_], sim_dict[cali]['y'][::thin_], s=2)
    plt.show()



## Specify the model parameters and values

parameter_files = [Path(os.path.join(cali_dir, 'params.yml')) for cali_dir, cali in zip(cali_dirs, cali_list)]
parameter_filenames = {filename.stem: filename for filename in parameter_files}
print(parameter_filenames)
parameter_filenames = {cali: filename for cali, filename in zip(cali_list, parameter_files)}

print(parameter_files)
print(parameter_filenames)

## Set this to what you ran the simulations with

model_specification, model_parameters = load_parameter_files(
    filenames=parameter_filenames,
    parameters=[
        "Parameters:BAL_f_accretion"
    ],
    log_parameters=[],
    parameter_printable_names=[
        r"$f_{\rm acc, BAL}$"
    ],
    #parameter_limits: Optional[List[List[float]]] = None, # for parameter limits
)

print(model_parameters)
print(model_specification)

modelvalues = {}

## Choose a min and max x value to limit the emulator to when training

# these are mass bins (make range more specific)

# Entropy Profiles
min_x_for_emulator = -1.7
max_x_for_emulator = 1.1

thin_ = 1
# thin_ = 10 # for CSHF

# rel_err = 0.02 # for Mbh-M*
rel_err = 0.05 # for CSFH

for cali in model_parameters.model_parameters:
    independent = sim_dict[cali]['x'][::thin_]
    dependent = sim_dict[cali]['y'][::thin_]
    
    # Must choose whether to specify own errorbars
    # or use maximum or mean errorbars (as swift emulator can only handle single value for errorbars)
#     dependent_error = np.abs(rel_err * dependent)
#    dependent_error = sim_dict[cali]['max_y_err'][::thin_]
    dependent_error = np.array([0]*len(dependent))

    range_condition = np.logical_and(independent>=min_x_for_emulator, independent<=max_x_for_emulator)
    
    modelvalues[cali] = {"independent": independent[range_condition],
                         "dependent": dependent[range_condition],
                         "dependent_error": dependent_error[range_condition]}

model_values = se.ModelValues(model_values=modelvalues)

## To check if there are any x bin values appearing only once
    # emulator needs more than one sample, cause it needs to compare
## This can make errors pop up when training the emulator, doing the parameter sensitivity analysis
## and performing the cross-check validation

comb_ = np.array([])
for key, val in model_values.model_values.items():
#     print(np.all(np.isfinite(val['independent'])))
#     print(np.all(np.isfinite(val['dependent'])))
#     print(np.all(np.isfinite(val['dependent_error'])))
#     print(max(val['independent']))
    comb_ = np.concatenate((comb_, val['independent']))
#     comb_ = np.append(comb_, val['independent'])
#     comb_ = np.append(comb_, max(val['independent']))
    
# Find unique values and their counts
unique_values, counts = np.unique(comb_, return_counts=True)
print("Unique x values:", unique_values)
print("Number of times:", counts)

# Find values that appear only once
values_appearing_once = unique_values[counts == 1]

# Print the values that appear only once



## Model Parameter Values

print(model_parameters)
print(model_specification)
model_parameters.plot_model(model_specification)


## TRAIN EMULATOR

for cali in model_parameters.model_parameters:
    independent = sim_dict[cali]['x'][::thin_]
    dependent = sim_dict[cali]['y'][::thin_]
    dependent_error = sim_dict[cali]['max_y_err'][::thin_]

    range_condition = np.logical_and(independent >= min_x_for_emulator,
                                     independent <= max_x_for_emulator)

    # Select values in range
    independent = independent[range_condition]
    dependent = dependent[range_condition]
    dependent_error = dependent_error[range_condition]

    # Remove NaNs and Infs
    mask = np.isfinite(independent) & np.isfinite(dependent) & np.isfinite(dependent_error)
    independent = independent[mask]
    dependent = dependent[mask]
    dependent_error = dependent_error[mask]

    # Save clean data
    modelvalues[cali] = {
        "independent": independent,
        "dependent": dependent,
        "dependent_error": dependent_error,
    }

emulator = gaussian_process.GaussianProcessEmulator()  # Default, MAKE IT

emulator.fit_model(model_specification=model_specification, #TRAIN IT
                            model_parameters=model_parameters,
                            model_values=model_values)
emulator.model_specification.sim_info = sim_info # LINK IT

#ALL TRAINING IS DONE HERE, ALL AFTER IS TO CHECK HOW WELL IT WORKED

## MAKE PREDICTIONS AND COMPARE:


for cali in emulator.model_parameters.model_parameters:
    pred_params = emulator.model_parameters[cali]
    pred_x = emulator.model_values[cali]['independent']
    sim_y = emulator.model_values[cali]['dependent']
    sim_yerr = emulator.model_values[cali]['dependent_error']

    pred_y, pred_var = emulator.predict_values(pred_x, pred_params)
    
    print("x: ", pred_x)
    print("y: ", pred_y)
    print("yerr: ", pred_var)
    print("Slope: ", [(y2-y1)/(x2-x1) for x1, x2, y1, y2 in zip(pred_x[:-1], pred_x[1:], pred_y[:-1], pred_y[1:])])

    plt.clf()
    plt.xlabel(sim_info['X_LABEL'])
    plt.ylabel(sim_info['Y_LABEL'])
    plt.title(r'%s, %s' % (cali, snap_prefix))
    plt.grid()
    plt.minorticks_on()
    
    plt.errorbar(pred_x, sim_y, color="black", ls=":", label="Simulation (Training Data)")
#     plt.scatter(pred_x, sim_y, color="black", label="Simulation (Training Data)")
    plt.fill_between(pred_x, sim_y-sim_yerr, sim_y+sim_yerr, color="black", alpha=0.2)
    
    plt.errorbar(pred_x, pred_y, yerr=pred_var, ls='-', label='Emulator')
#     plt.scatter(pred_x, pred_y, label='Emulator')
    plt.legend()
    plt.show()

## See how far off each prediction is from its simulation

pred_over_sim_list = []

for cali in emulator.model_parameters.model_parameters:
    pred_params = emulator.model_parameters[cali]
    pred_x = emulator.model_values[cali]['independent']
    sim_y = emulator.model_values[cali]['dependent']

    pred_y, pred_var = emulator.predict_values(pred_x, pred_params)
    
    pred_over_sim = pred_y/sim_y
    pred_over_sim_list = np.append(pred_over_sim_list, pred_over_sim)

    plt.plot(pred_x, pred_over_sim)
#     plt.scatter(pred_x, pred_over_sim)

plt.grid()
plt.minorticks_on()
plt.xlabel(sim_info['X_LABEL'])
plt.ylabel("Prediction / Simulation")
plt.show()


## PERFORM CROSS CHECK (VALIDATION)

emulator_ccheck = cross_check.CrossCheck(hide_progress=False)
emulator_ccheck.build_emulators(model_specification=model_specification,
                        model_parameters=model_parameters,
                        model_values=model_values)

data_by_cc = emulator_ccheck.build_mocked_model_values_original_independent()


cc_list = []
cc_over_og_list = []

for unique_identifier in emulator_ccheck.model_values.model_values:
    cc_over_og = data_by_cc[unique_identifier]["dependent"] / \
                model_values[unique_identifier]["dependent"]

    cc_list = np.append(cc_list, data_by_cc[unique_identifier]["dependent"])
    cc_over_og_list = np.append(cc_over_og_list, cc_over_og)
    
    plt.plot(data_by_cc[unique_identifier]["independent"], cc_over_og)

plt.grid()
plt.minorticks_on()
plt.xlabel(sim_info['X_LABEL'])
plt.ylabel("Cross-check / Truth")
plt.show()

## Cross-check does better if y-axis value is negative (strange, but good. Machine learning is working....?)
## Original prediction does better if y-axis value is positive

pred_list = []

for unique_identifier in emulator_ccheck.model_values.model_values:
    
    pred_params = emulator.model_parameters[unique_identifier]
    pred_x = emulator.model_values[unique_identifier]['independent']
    sim_y = emulator.model_values[unique_identifier]['dependent']

    pred_y, pred_var = emulator.predict_values(pred_x, pred_params)
    
    pred_list = np.append(pred_list, pred_y)
    
    
    val_ = np.log10(np.abs((data_by_cc[unique_identifier]["dependent"] - sim_y)/(pred_y - sim_y)))
    plt.plot(data_by_cc[unique_identifier]["independent"], val_)

plt.grid()
plt.minorticks_on()
plt.xlabel(sim_info['X_LABEL'])
plt.ylabel(r'$\log{ \left| \frac{Cross \, check - Simulation}{Prediction - Simulation} \right|}$')
plt.show()


## Second way of comparing the cross-check for each simulation to the prediction as found by
## the fully trained emulator

plt.plot([min(pred_list), max(pred_list)], [min(pred_list), max(pred_list)], color='black')
plt.scatter(pred_list, cc_list, s=1)

plt.grid()
plt.minorticks_on()
plt.xlabel("Prediction")
plt.ylabel("Cross-check")
plt.show()


emulator_ccheck.plot_results(emulate_at=data_by_cc[unique_identifier]["independent"],
                         xlabel=sim_info['X_LABEL'],
                         ylabel=sim_info['Y_LABEL'])

total_mean_squared = emulator_ccheck.get_mean_squared()[0]
print("Total mean squared of entire set of left-out simulations: %s" % total_mean_squared)

## As in FLAMINGO paper (Kugel+22), error on emulator approximated as standard deviation of
## cross-check values divided by original simulation values

# sigma_ccheck = np.nanstd(cc_over_og_list)#, ddof=0)
# sigma_ccheck = total_mean_squared
sigma_ccheck = 0.1 #set to 10 percent error, easiest lol
print('sigma_ccheck =', sigma_ccheck)
emulator.model_specification.sigma_ccheck = sigma_ccheck

def emulator_model(x, pred_params, emulator):
    pred_y, pred_var = emulator.predict_values(x, pred_params)
    
    pred_std = np.sqrt(pred_var)
    pred_ccheck_std = np.abs(emulator.model_specification.sigma_ccheck * pred_y)
    pred_total = np.sqrt(pred_std**2 + pred_ccheck_std**2)
    
    return pred_y, pred_ccheck_std


# Make predictions with trained emulator (but now with emulator uncertainties!) to compare to input

for cali in emulator.model_parameters.model_parameters:
    pred_params = emulator.model_parameters[cali]
    pred_x = emulator.model_values[cali]['independent']
    sim_y = emulator.model_values[cali]['dependent']
    sim_yerr = emulator.model_values[cali]['dependent_error']
    
    pred_y, pred_err = emulator_model(pred_x, pred_params, emulator)

    print("x: ", pred_x)
    print("pred_y: ", pred_y)
    print("pred_yerr: ", pred_err) 

    plt.clf()
    plt.xlabel(sim_info['X_LABEL'])
    plt.ylabel(sim_info['Y_LABEL'])
    plt.title(r'%s, %s' % (cali, snap_prefix))
    plt.grid()
    plt.minorticks_on()
    
#     plt.errorbar(pred_x, sim_y, yerr=sim_yerr, color="black", ls=":", label="Simulation (Training Data)")
    plt.errorbar(pred_x, sim_y, color="black", ls=":", label="Simulation (Training Data)")
#     plt.scatter(pred_x, sim_y, color="black", label="Simulation (Training Data)")
    plt.fill_between(pred_x, sim_y-sim_yerr, sim_y+sim_yerr, color='black', alpha=0.2)

#     plt.errorbar(pred_x, pred_y, yerr=pred_err, ls='-', fmt='o', label='Emulator')
    plt.errorbar(pred_x, pred_y, yerr=None, ls='-', label='Emulator')
#     plt.scatter(pred_x, pred_y, label='Emulator')
    plt.fill_between(pred_x, pred_y-pred_err, pred_y+pred_err, alpha=0.3)
    plt.legend()
    plt.show()



## SWEEPS OF PARAMETER SPACE

# Get full range of x values that have been used in training emulator

#another way of checking parameter sensitivity

sim_x = []
for cali in emulator.model_values.model_values:
    sim_x_curr = emulator.model_values[cali]['independent']
    for x_val in sim_x_curr:
        if x_val not in sim_x:
            sim_x.append(x_val)
            
sim_x = np.array(np.sort(sim_x))
print(sim_x)

# Arbitrarily choose initial parameters to sweep from
# Doesn't matter too much though, because it ends up covering
# roughly the entire range of parameters used

centre = emulator.model_parameters['run01']

for ii in range(len(emulator.parameter_order)):
    param_name = emulator.parameter_order[ii]
    param_printable_name = emulator.model_specification.parameter_printable_names[ii]
    
    Mock_values, Mock_parameters = mock_sweep(emulator, emulator.model_specification, 6, param_name, centre)

    plt.clf()
    
    for mock_name in Mock_values.keys():
        plt.plot(Mock_values[mock_name]["independent"],
                Mock_values[mock_name]["dependent"],
                label = "%s = %.4g" % (param_printable_name, Mock_parameters[mock_name][param_name]))
#         plt.scatter(Mock_values[mock_name]["independent"],
#                 Mock_values[mock_name]["dependent"],
#                 label = "%s = %.4g" % (param_printable_name, Mock_parameters[mock_name][param_name]))
    
    plt.xlabel(sim_info['X_LABEL'])
    plt.ylabel(sim_info['Y_LABEL'])
    plt.grid()
    plt.minorticks_on()
    plt.legend()
    plt.show()




## CHECKING HYPERPARAMETERS


emulator.kernel.get_parameter_dict(include_frozen=True)
emulator_mcmc = gaussian_process_mcmc.GaussianProcessEmulatorMCMC(burn_in_steps=1, mcmc_steps=1000)#, hide_progress=False)
emulator_mcmc.fit_model(model_specification=model_specification,
                        model_parameters=model_parameters,
                        model_values=model_values)

emulator_mcmc.kernel.get_parameter_dict(include_frozen=True)


# save emulator

def save_object(obj, filename):
    with open(filename, 'wb') as f:  # Overwrites any existing file.
        dill.dump(obj, f, dill.HIGHEST_PROTOCOL)

emulator_filename = 'emulator_'+obs_name+'.pkl'
emulator_path = os.path.join(root_dir, 'emulators')
if not os.path.isdir(emulator_path):
    os.makedirs(emulator_path, mode = 0o755, exist_ok=True)
save_object(emulator, os.path.join(emulator_path, emulator_filename))