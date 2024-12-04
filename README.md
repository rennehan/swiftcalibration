# swiftcalibration
Calibrate Swift

# Instructions for use

## Creating calibration simulations
- Create a python environment using one of the requirements.txt files. Activate it.

- Copy required simulation files to the "data" directory (you can remove whatever is not necessary for you already in there). Make sure to write the names of these files in the files_to_link list in "template/generate_calibrations.py".

- Alter "template/design.py" to include desired parameters to calibrate, their ranges, their printable names, whether they should be sampled in log space, how many simulations you would like to run (32 tends to be good), and the base parameter file you are using (which must be compatible with the swift executable you are using)

- cd to "template", run "python design.py" and then "python generate_calibrations.py <system_name> </path/to/initial_conditions_file>".

- Run "python submit_calibrations_auto.py" if you want the simulations to restart if they fail; if not, run "python submit_calibrations.py" first and then "python submit_calibrations_restart.py" subsequent times after they fail or hit their wall time.

## Performing Calibration
- There are scripts for generating caesar files and some basic observables in "scripts", which can all be run together by submitting the job script "job_gen_observables.sh". This will produce velociraptor hdf5 files for the observables in each calibration's directory.

- The jupyter notebook "gen_obs_data.ipynb" can be used to produce velociraptor hdf5 files of observational data, which can be stored wherever desired.

- You can then open the jupyter notebook "gen_swift_emulator.ipynb", and run through all the cells to generate and save a different emulator for each observable.

- Finally, the jupyter notebook "swift_emulator_joint_mcmc.ipynb" is used to find the best of the calibration simulations, and then use the emulators for each observable jointly in an MCMC to find the overall best-fit parameters.

- If desired, a single full simulation can be run with these best-fit parameters, and then scripts/gen_sim_data.ipynb can be used to look at the observables of that simulation.
