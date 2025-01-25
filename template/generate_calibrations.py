import os
import numpy as np
from glob import glob
from pathlib import Path
import argparse as ap

#import logging
#logger = logging.getLogger(__name__)

def find_git_dir(directory):
    "Find the correct git dir; move upwards if .git folder is not found here"
    absdir = os.path.abspath(directory)
    gitdir = os.path.join(absdir, ".git")
    if os.path.isdir(gitdir):
        return gitdir
    parentdir = os.path.dirname(absdir)
    if absdir == parentdir:
        # We reached root and found no gitdir
#        logger.warning("No git dir found")
        print('No git dir found')

        return None
    return find_git_dir(parentdir)

git_top_dir = find_git_dir('.')[:-5]


# ics_file = The HDF5 file with the ICs.
parser = ap.ArgumentParser()
parser.add_argument("ics_file")
args = parser.parse_args()


data_dir = '../../../data'
print(data_dir)

files_to_link = ["yieldtables",
                 "photometry",
                 "output_list_cali.txt",
                 "chem5",
                 "coolingtables",
                 "CloudyData_UVB=FG2011_shielded.h5"]#,
#                 args.ics_file]



parameter_files = [Path(x) for x in glob("./original_ymls/*.yml")]
parameter_filenames = {filename.stem: filename for filename in parameter_files}

cali_dir = './calibrations'
if not os.path.isdir(cali_dir):
    os.makedirs(cali_dir, mode = 0o755, exist_ok = True)

new_ymls = []
for v,k in enumerate(parameter_filenames):
    cali_name = 'cali_%04d' % int(k)
    cali_path = os.path.join(cali_dir, cali_name)
    yml_file = 'cali_%04d.yml' % int(k)

    if os.path.isdir(cali_path):
        os.system('rm -rf %s' % cali_path)
    os.makedirs(cali_path, mode = 0o755, exist_ok = True)
    os.makedirs('%s/src' % cali_path, mode = 0o755, exist_ok = True)
    os.system('rsync -av ~/src/swiftsim/* %s/src/' % cali_path)
    os.makedirs('%s/src/.git' % cali_path, mode = 0o755, exist_ok = True)
    os.system('rsync -av ~/src/swiftsim/.git/* %s/src/.git/' % cali_path)
    os.system('ln -s `pwd`/%s/src/swift %s/swift' % (cali_path, cali_path))

    os.system('ln -s %s/%s %s/%s' % (data_dir, args.ics_file, cali_path, 'ics_file.hdf5'))
    for file_to_link in files_to_link:
        os.system('ln -s %s/%s %s/%s' % (data_dir, file_to_link, cali_path, file_to_link))

    os.system('cp ./original_ymls/%d.yml %s' % (int(k), os.path.join(cali_path, yml_file)))


os.system('python ./generate_jobs.py')
