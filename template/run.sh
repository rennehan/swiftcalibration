## Run this master script, i.e. <bash run.sh> or <./run.sh>
## You must have the correct python environment activated to run this script
## and set up the correct variables, paths, parameters, etc... in design_calibrations.py, generate_calibrations.py, move_ymls.py, generate_jobs.py, and submit_calibrations.py

source /scratch/aspadawe/calibration/pyenvs/swift-calibration-2/bin/activate

python design_calibrations.py
python generate_calibrations.py
# python move_ymls.py
python generate_jobs.py
# python submit_calibrations.py

