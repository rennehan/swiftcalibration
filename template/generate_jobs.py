import os
import numpy as np
from glob import glob
from pathlib import Path

#data_dir = '/home/b/babul/rennehan/work/swift'
data_dir = '/scratch/b/babul/aspadawe/swift_tests/s50n512_simba_ps2020_sphenix_calibration7_mpi'
#files_to_link = ['swift', 'yieldtables', 'photometry', 'output_list_cali.txt', 'UV_dust1_CR1_G1_shield1.hdf5', 'Jenny_N128L12.dat.hdf5', 'velo.cfg']
#files_to_link = ['swift', 'yieldtables', 'photometry', 'output_list_cali.txt', 'UV_dust1_CR1_G1_shield1.hdf5', 'Jenny_N128L12.dat.hdf5', 'chem5_tables', 'CloudyData_UVB=FG2011_shielded.h5']
files_to_link = ['swift_mpi_ps2020_sphenix_v2', 'yieldtables', 'photometry', 'output_list_cali.txt', 'UV_dust1_CR1_G1_shield1.hdf5', 'Jenny_N512L50.hdf5', 'chem5_tables', 'CloudyData_UVB=HM2012_shielded.h5']

with open('job.sh', 'r') as f:
    job_data = f.read()

with open('job_restart.sh', 'r') as f:
    job_restart_data = f.read()

with open('job_auto.sh', 'r') as f:
    job_auto_data = f.read()

parameter_files = [Path(x) for x in glob("./original_ymls/*.yml")]
parameter_filenames = {filename.stem: filename for filename in parameter_files}

cali_dir = './calibrations'
if not os.path.isdir(cali_dir):
    os.makedirs(cali_dir, mode = 0o755, exist_ok = True)

new_jobs = []
new_jobs_restart = []
new_jobs_auto = []
new_ymls = []
for v,k in enumerate(parameter_filenames):
    cali_name = 'cali_%04d' % int(k)
    cali_path = os.path.join(cali_dir, cali_name)
    yml_file = 'cali_%04d.yml' % int(k)
    job_file = 'job_%04d.sh' % int(k)
    job_restart_file = 'job_restart_%04d.sh' % int(k)
    job_auto_file = 'job_auto_%04d.sh' % int(k)


    new_jobs.append(job_data.replace('JOB_NAME', cali_name))
    new_jobs[-1] = new_jobs[-1].replace('YML_FILE', yml_file)
    new_jobs[-1] = new_jobs[-1].replace('JOB_FILE', job_file)
    with open(os.path.join(cali_path, job_file), 'w') as f:
        f.write(new_jobs[-1])

    new_jobs_restart.append(job_restart_data.replace('JOB_NAME', cali_name))
    new_jobs_restart[-1] = new_jobs_restart[-1].replace('YML_FILE', yml_file)
    new_jobs_restart[-1] = new_jobs_restart[-1].replace('JOB_FILE', job_restart_file)
    with open(os.path.join(cali_path, job_restart_file), 'w') as f:
        f.write(new_jobs_restart[-1])

    new_jobs_auto.append(job_auto_data.replace('JOB_NAME', cali_name))
    new_jobs_auto[-1] = new_jobs_auto[-1].replace('YML_FILE', yml_file)
    new_jobs_auto[-1] = new_jobs_auto[-1].replace('JOB_FILE', job_auto_file)
    with open(os.path.join(cali_path, job_auto_file), 'w') as f:
        f.write(new_jobs_auto[-1])
