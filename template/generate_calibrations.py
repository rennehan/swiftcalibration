import os
import numpy as np
from glob import glob
from pathlib import Path
import argparse

#import logging
#logger = logging.getLogger(__name__)


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



output_dir = '/scratch/aspadawe/snapshots/HyenasC/L1/SimbaC_L1_Calibration/halo_2205-correct_jet-calibrations/128_sims-6_params/'
cali_dir = '/scratch/aspadawe/sims/HyenasC/L1/SimbaC_L1_Calibration/trillium/gizmo-hyenasc-l1-correct_jet-vary_params/cali_params/128_sims-6_params/'
param_type = 'tex'
ics_file = '/project/rrg-babul-ad/wcui/HYENAS/ICs/level1/halo_2205'

# def find_git_dir(directory):
#     "Find the correct git dir; move upwards if .git folder is not found here"
#     absdir = os.path.abspath(directory)
#     gitdir = os.path.join(absdir, ".git")
#     if os.path.isdir(gitdir):
#         return gitdir
#     parentdir = os.path.dirname(absdir)
#     if absdir == parentdir:
#         # We reached root and found no gitdir
# #        logger.warning("No git dir found")
#         print('No git dir found')

#         return None
#     return find_git_dir(parentdir)

# git_top_dir = find_git_dir('.')[:-5]


# ics_file = The HDF5 file with the ICs.
# parser = argparse.ArgumentParser()
# parser.add_argument("ics_file")
# args = parser.parse_args()


# data_dir = '../../../data'
# print(data_dir)

files_to_link = {
    'GIZMO_EXE': '/scratch/aspadawe/sims/HyenasC/L1/SimbaC_L1_Calibration/trillium/gizmo-hyenasc-l1-correct_jet-vary_params/GIZMO_EXE',
    'ics_file': ics_file,
    'param_files': '/scratch/aspadawe/sims/HyenasC/L1/SimbaC_L1_Calibration/trillium/gizmo-hyenasc-l1-correct_jet-vary_params/param_files/',
}

# required_files = [Path(x) for x in glob(f'/scratch/aspadawe/snapshots/required_files/*')]
required_files = [x for x in glob(f'/scratch/aspadawe/snapshots/required_files/*')]
for file in required_files:
    # files_to_link[file.stem] = file
    files_to_link[file.split('/')[-1]] = file

print()
print(files_to_link)
print()

# files_to_link = ['/scratch/aspadawe/sims/HyenasC/L1/SimbaC_L1_Calibration/trillium/gizmo-hyenasc-l1-correct_jet-vary_params/GIZMO_EXE',
#                  '/scratch/aspadawe/sims/HyenasC/L1/SimbaC_L1_Calibration/trillium/gizmo-hyenasc-l1-correct_jet-vary_params/param_files/',
#                  '/scratch/aspadawe/snapshots/required_files/*',]
# files_to_link = ["yieldtables",
#                  "photometry",
#                  "output_list_cali.txt",
#                  "chem5",
#                  "coolingtables",
#                  "CloudyData_UVB=FG2011_shielded.h5",
#                  "snapshot_BAL_0001.hdf5",
#                  "snapshot_BAL_0002.hdf5",
#                  "snapshot_BAL_0003.hdf5",
#                  "snapshot_BAL_0004.hdf5",
#                  "snapshot_BAL_0005.hdf5",
#                  "snapshot_BAL_0006.hdf5",
#                  "snapshot_BAL_0007.hdf5",
#                  "snapshot_BAL_0008.hdf5",
#                  "snapshot_BAL_0009.hdf5",
#                  "snapshot_BAL_0010.hdf5",
#                  ]#,
#                 args.ics_file]



# parameter_files = [Path(x) for x in glob("./original_ymls/*.yml")]
parameter_files = [Path(x) for x in glob(f"{cali_dir}/*.{param_type}")]
parameter_filenames = {filename.stem: filename for filename in parameter_files}

# cali_dir = './calibrations'
if not os.path.isdir(output_dir):
    os.makedirs(output_dir, mode = 0o755, exist_ok = True)

# new_paramfiles = []
for v,k in enumerate(parameter_filenames):
    # cali_name = 'cali_%04d' % int(k)
    cali_name = f'cali_{int(k):04d}'
    cali_path = os.path.join(output_dir, cali_name)
    # param_file = 'cali_%04d.%s' % (int(k), args.param_type)
    param_file = f'cali_{int(k):04d}.{param_type}'

    if os.path.isdir(cali_path):
        # os.system('rm -rf %s' % cali_path)
        os.system(f'rm -rf {cali_path}')
    os.makedirs(cali_path, mode = 0o755, exist_ok = True)
    os.makedirs(os.path.join(cali_path, 'slurm_files'), mode = 0o755, exist_ok = True)

    # os.system(f'ln -s ')

    # os.makedirs('%s/src' % cali_path, mode = 0o755, exist_ok = True)
    # os.system('rsync -av ~/src/swiftsim/* %s/src/' % cali_path)
    # os.makedirs('%s/src/.git' % cali_path, mode = 0o755, exist_ok = True)
    # os.system('rsync -av ~/src/swiftsim/.git/* %s/src/.git/' % cali_path)
    # os.system('ln -s `pwd`/%s/src/swift %s/swift' % (cali_path, cali_path))

    # os.system('ln -s %s/%s %s/%s' % (data_dir, args.ics_file, cali_path, 'ics_file.hdf5'))

    for file_name, file_to_link in files_to_link.items():
        os.system(f'ln -s {file_to_link} {os.path.join(cali_path, file_name)}')

    # os.system('cp -rf ./original_ymls/%d.yml %s' % (int(k), os.path.join(cali_path, yml_file)))
    os.system(f'cp -rf {cali_dir}/{int(k):d}.{param_type} {os.path.join(cali_path, param_file)}')


# os.system('python ./generate_jobs.py')
