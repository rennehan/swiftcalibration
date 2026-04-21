import os
import numpy as np
from glob import glob
from pathlib import Path
import argparse

#import logging
#logger = logging.getLogger(__name__)


print("\nGenerating calibrations...\n")


## Command line keyword arguments
# parser = argparse.ArgumentParser(description="Generate calibrations.")
# # parser.add_argument('--code_dir', action='store', type=str, required=True, 
# #                    help='Path to directory containing executable gizmo file.')
# # parser.add_argument('--output_dir', action='store', type=str, required=True, 
# #                    help='Path to directory where calibration output will be written.')
# # parser.add_argument('--cali_dir', action='store', type=str, required=True, 
# #                    help='Path to directory where calibration parameter files are stored.')
# # parser.add_argument('--param_type', action='store', type=str, default='yml', choices=['yml','tex'],
# #                    help='Type of parameter file to use for calibration (yml or tex).')
# # parser.add_argument('--ics_file', action='store', type=str, required=True, 
# #                    help='Path to initial conditions file (HDF5 format).')
# # args = parser.parse_args()



## Specify paths and parameters for generating calibrations; modify as needed for your calibration parameters and model
output_dir = '/path/to/output/calibration/directory/'  # directory where calibration output subdirectories will be written
cali_dir = '/path/to/calibration/parameter/files/'  # same as output_dir in design_calibrations.py
param_type = 'tex'  # type of parameter file to use for calibration (yml or tex)
model_info_filename = 'model_info.pkl'  # same as model_info_filename in design_calibrations.py


## Define files to link into each calibration subdirectory
## Currently for gizmo-simba-hyenasc; modify as needed for your calibration parameters and model
ics_file = '/path/to/initial/conditions/file'
files_to_link = {
    'GIZMO_EXE':'/path/to/gizmo/executable',
    'ics_file':ics_file,
    'param_files':'/scratch/aspadawe/sims/HyenasC/L1/SimbaC_L1_Calibration/trillium/gizmo-hyenasc-l1-correct_jet-vary_params/param_files/',
}

# required_files = [Path(x) for x in glob(f'/scratch/aspadawe/snapshots/required_files/*')]
required_files = [x for x in glob(f'/scratch/aspadawe/snapshots/required_files/*')]
for file in required_files:
    # files_to_link[file.stem] = file
    files_to_link[file.split('/')[-1]] = file

print("\tWill link the following files into each calibration subdirectory:")
print(files_to_link)


## Check if output directory exists; if not, create it
if not os.path.isdir(output_dir):
    print(f"\tOutput directory {output_dir} does not exist; creating it.")
    os.makedirs(output_dir, mode = 0o755, exist_ok = True)


## Copy model_info file to output directory
print(f"\tCopying model information file {model_info_filename} to output directory {output_dir}...")
os.system(f'cp -rf {os.path.join(cali_dir, model_info_filename)} {os.path.join(output_dir, model_info_filename)}')


## Find parameter files in cali_dir with specified type (yml or tex) and generate calibration setups for each parameter file
# parameter_files = [Path(x) for x in glob("./original_ymls/*.yml")]
# parameter_files = [Path(x) for x in glob(f"{cali_dir}/*.{param_type}")]
# parameter_files = [Path(x) for x in glob(os.path.join(cali_dir, f"cali_*.{param_type}"))]  # look for parameter files in cali_dir with format cali_*.yml or cali_*.tex
parameter_files = [Path(x) for x in glob(os.path.join(cali_dir, f"*.{param_type}"))]  # look for parameter files in cali_dir with format *.yml or *.tex
parameter_filenames = {filename.stem:filename for filename in parameter_files}
print(f"\tFound {len(parameter_files)} parameter files in {cali_dir} with type {param_type}")

cali_dict = {
    cali:{
        'cali_dir':cali_path,
        'param_file':f'params.{param_type}',
    } 
    for cali, cali_path in parameter_filenames.items()
}

# new_paramfiles = []
print(f"\tGenerating calibration setups for each parameter file...")
for i, (cali_num, cali_input_path) in enumerate(parameter_filenames.items()):
    cali_name = f'cali_{int(cali_num):04d}'
    cali_output_path = os.path.join(output_dir, cali_name)

    cali_yml_input_path = os.path.join(cali_dir, f'{cali_num}.yml')

    ## Set parameter file name for calibration; e.g. cali_0001.tex, cali_0001.yml,
    ## or just params.tex or params.yml if you want the same parameter file name for each calibration
    # param_file = f'{cali_name}.{param_type}'
    # param_file = f'params.{param_type}'
    param_file = 'params'

    if os.path.isdir(cali_output_path):
        print(f"\tCalibration output directory {cali_output_path} already exists; removing this calibration to create a new one.")
        os.system(f'rm -rf {cali_output_path}')
    os.makedirs(cali_output_path, mode = 0o755, exist_ok = True)

    ## Make subdirectory for output files, e.g. slurm files, if your setup uses slurm
    os.makedirs(os.path.join(cali_output_path, 'slurm_files'), mode = 0o755, exist_ok = True)

    ## Link required files into calibration subdirectory
    # print(f"\tLinking required files into calibration subdirectory {cali_output_path}...")
    for file_name, file_to_link in files_to_link.items():
        os.system(f'ln -s {file_to_link} {os.path.join(cali_output_path, file_name)}')

    ## Copy parameter file to calibration subdirectory
    # os.system(f'cp -rf {cali_dir}/{int(cali_num):d}.{param_type} {os.path.join(cali_output_path, param_file)}')
    os.system(f'cp -rf {cali_input_path} {os.path.join(cali_output_path, f"{param_file}.{param_type}")}')

    ## If parameter file is not already in .yml format, copy the .yml version of the parameter file to the calibration subdirectory as well,
    ## since this is what the swiftemulator functions need
    if param_type.lower() != 'yml':
        os.system(f'cp -rf {cali_yml_input_path} {os.path.join(cali_output_path, f"{param_file}.yml")}')

    print(f"\tGenerated calibration setup {i+1}/{len(parameter_filenames)}: {cali_output_path}")


print("\nFinished generating calibrations.\n")