import os
import numpy as np
from glob import glob
from pathlib import Path

#data_dir = '/home/b/babul/rennehan/work/swift'
data_dir = '/scratch/b/babul/aspadawe/swift_tests/s25n256_simba_cloudy'
#files_to_link = ['swift', 'yieldtables', 'photometry', 'output_list_cali.txt', 'UV_dust1_CR1_G1_shield1.hdf5', 'Jenny_N128L12.dat.hdf5', 'velo.cfg']
#files_to_link = ['swift', 'yieldtables', 'photometry', 'output_list_cali.txt', 'UV_dust1_CR1_G1_shield1.hdf5', 'Jenny_N128L12.dat.hdf5', 'chem5_tables', 'CloudyData_UVB=FG2011_shielded.h5']
files_to_link = ['swift', 'yieldtables', 'photometry', 'output_list_cali.txt', 'UV_dust1_CR1_G1_shield1.hdf5', 'Jenny_N256L25.dat.hdf5', 'chem5_tables', 'CloudyData_UVB=HM2012_shielded.h5']

with open('job.sh', 'r') as f:
    job_data = f.read()

with open('job_restart.sh', 'r') as f:
    job_restart_data = f.read()

parameter_files = [Path(x) for x in glob("./original_ymls/*.yml")]
parameter_filenames = {filename.stem: filename for filename in parameter_files}

cali_dir = './calibrations'
if not os.path.isdir(cali_dir):
    os.makedirs(cali_dir, mode = 0o755, exist_ok = True)

new_jobs = []
new_jobs_restart = []
new_ymls = []
for v,k in enumerate(parameter_filenames):
    cali_name = 'cali_%04d' % int(k)
    cali_path = os.path.join(cali_dir, cali_name)
#    yml_file = 'swimba_%04d.yml' % int(k)
    yml_file = 'cali_%04d.yml' % int(k)
    job_file = 'job_%04d.sh' % int(k)
    job_restart_file = 'job_restart_%04d.sh' % int(k)

    if os.path.isdir(cali_path):
        os.system('rm -rf %s' % cali_path)
    os.makedirs(cali_path, mode = 0o755, exist_ok = True)

    for file_to_link in files_to_link:
        os.system('ln -s %s/%s %s/%s' % (data_dir, file_to_link, cali_path, file_to_link))

    new_jobs.append(job_data.replace('JOB_NAME', cali_name))
    new_jobs[-1] = new_jobs[-1].replace('YML_FILE', yml_file)
    with open(os.path.join(cali_path, job_file), 'w') as f:
        f.write(new_jobs[-1])

    new_jobs_restart.append(job_restart_data.replace('JOB_NAME', cali_name))
    new_jobs_restart[-1] = new_jobs_restart[-1].replace('YML_FILE', yml_file)
    with open(os.path.join(cali_path, job_restart_file), 'w') as f:
        f.write(new_jobs_restart[-1])

    os.system('cp ./original_ymls/%d.yml %s' % (int(k), os.path.join(cali_path, yml_file)))
