# swiftcalibration
Calibrate Swift-based (and GADGET/GIZMO-based) simulations.

# Instructions for use

First, clone this repository.

## Creating calibration simulations
- Create a python environment using one of the requirements.txt files. Activate it.

- Copy required simulation files to the "data" directory (you can remove whatever is not necessary for you already in there). Make sure to write the names of these files in the files_to_link list in "template/generate_calibrations.py".

- Alter "template/design.py": Set the number of calibration simulations to run (<num_simulations>), the path to the base/template parameter file (<base_param_file>), and the path to the directory where the calibration parameter files will be written (<output_path>). The code is designed for SWIFT-based simulations, but a work-around has been implemented for GADGET/GIZMO-based simulations; if the latter is the case, set both <convert_input_tex_to_yml> and <convert_output_yml_to_tex> to True. Additionally, set the (number of) desired parameters to calibrate, their ranges, their printable names, and whether they should be sampled in log space.

- Alter "template/generate_calibrations.py": Set the path to the directory where the simulation outputs will be written (<output_dir>), the path to the directory containing the calibration parameter files (<cali_dir>, this will be the same as <output_path> in "template/design.py"), the extension of the parameter file (<param_type>, if the conversion between yml and tex/txt param files was used in design.py, this is important), and the path to the initial conditions (ICs) file for the simulation (<ics_file>). Also add the paths for all files required to run the simulation to the files_to_link dictionary, using the dictionary keys to specify what these files are called in the parameters file.

- Alter "template/generate_jobs.py": 

- cd to "template", run "python design.py" and then "python generate_calibrations.py <system_name> </path/to/initial_conditions_file>".

- Run "python submit_calibrations_auto.py" if you want the simulations to restart if they fail; if not, run "python submit_calibrations.py" first and then "python submit_calibrations_restart.py" subsequent times after they fail or hit their wall time.

## Performing Calibration
- There are scripts for generating caesar files and some basic observables in "scripts", which can all be run together by submitting the job script "job_gen_observables.sh". This will produce velociraptor hdf5 files for the observables in each calibration's directory.

- The jupyter notebook "gen_obs_data.ipynb" can be used to produce velociraptor hdf5 files of observational data, which can be stored wherever desired.

- You can then open the jupyter notebook "gen_swift_emulator.ipynb", and run through all the cells to generate and save a different emulator for each observable.

- Finally, the jupyter notebook "swift_emulator_joint_mcmc.ipynb" is used to find the best of the calibration simulations, and then use the emulators for each observable jointly in an MCMC to find the overall best-fit parameters.

- If desired, a single full simulation can be run with these best-fit parameters, and then scripts/gen_sim_data.ipynb can be used to look at the observables of that simulation.
