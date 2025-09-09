# swiftcalibration
Modified to calibrate GIZMO-based simulations, specifically for entropy data.

# Instructions for use

First, clone this repository.

## Step-by-Step Walkthrough

- For GIZMO simulations, run the parameter file from the simulation through "paramfile\_tex\_to\_yml.py", which will convert it from .tex to .yml format. Within this script, change the input\_tex\_filepath and output\_yml\_filepath to your desired filepaths.
  
- Run the "gen\_entr.py" script, with your desired filepaths inserted. The script will loop through all your simulated runs and generate entropy profile .hdf5 files, which will be used to train the emulator.
  
- You can then open the jupyter notebook "gen\_swift\_emulator.ipynb", and run through all the cells to generate and save a different emulator for each observable.
  
- Now, run "gen\_entr\_observational.py" to create an .hdf5 file for observational entropy profiles. You will need to input a mass range and an observational catalog (ACCEPT, CLoGS, or Sun+2009) to take data from.
    ~ If your simulation halo is ~M500 = 13, then it would be logical to pick an observational mass range of about 12.8-13.2.
  
- Use the observational entropy profile file you've created along with the trained emulator (.pkl file) and your snapshot information to run through the run through the jupyter notebook "swift\_emulator\_joint\_mcmc.ipynb". This notebook is used to find the best of the calibration simulations, and then use the emulators for each observable jointly in an MCMC to find the overall best-fit parameters. 

Email spencerlockwood@uvic.ca for any questions
