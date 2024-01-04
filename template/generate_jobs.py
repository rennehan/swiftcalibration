import os
import numpy as np
from glob import glob
from pathlib import Path
import argparse as ap

# System = niagara, rusty
parser = ap.ArgumentParser()
parser.add_argument("system")
args = parser.parse_args()


with open('job.sh', 'r') as f:
    job_data = f.read()

with open('job_restart.sh', 'r') as f:
    job_restart_data = f.read()

with open('job_mpi.sh', 'r') as f:
    job_mpi_data = f.read()

with open('job_restart_mpi.sh', 'r') as f:
    job_restart_mpi_data = f.read()

with open('job_auto.sh', 'r') as f:
    job_auto_data = f.read()

parameter_files = [Path(x) for x in glob("./original_ymls/*.yml")]
parameter_filenames = {filename.stem: filename for filename in parameter_files}

cali_dir = './calibrations'
if not os.path.isdir(cali_dir):
    os.makedirs(cali_dir, mode = 0o755, exist_ok = True)

new_jobs = []
new_jobs_restart = []
new_jobs_mpi = []
new_jobs_restart_mpi = []
new_jobs_auto = []
new_ymls = []
for v,k in enumerate(parameter_filenames):
    cali_name = 'cali_%04d' % int(k)
    cali_path = os.path.join(cali_dir, cali_name)
    yml_file = 'cali_%04d.yml' % int(k)
    job_file = 'job_%04d.sh' % int(k)
    job_restart_file = 'job_restart_%04d.sh' % int(k)
    job_mpi_file = 'job_mpi_%04d.sh' % int(k)
    job_restart_mpi_file = 'job_restart_mpi_%04d.sh' % int(k)
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

    new_jobs_mpi.append(job_mpi_data.replace('JOB_NAME', cali_name))
    new_jobs_mpi[-1] = new_jobs_mpi[-1].replace('YML_FILE', yml_file)
    new_jobs_mpi[-1] = new_jobs_mpi[-1].replace('JOB_FILE', job_mpi_file)
    with open(os.path.join(cali_path, job_mpi_file), 'w') as f:
        f.write(new_jobs_mpi[-1])

    new_jobs_restart_mpi.append(job_restart_mpi_data.replace('JOB_NAME', cali_name))
    new_jobs_restart_mpi[-1] = new_jobs_restart_mpi[-1].replace('YML_FILE', yml_file)
    new_jobs_restart_mpi[-1] = new_jobs_restart_mpi[-1].replace('JOB_FILE', job_restart_mpi_file)
    with open(os.path.join(cali_path, job_restart_mpi_file), 'w') as f:
        f.write(new_jobs_restart_mpi[-1])
    
    new_jobs_auto.append(job_auto_data.replace('JOB_NAME', cali_name))
    new_jobs_auto[-1] = new_jobs_auto[-1].replace('YML_FILE', yml_file)
    new_jobs_auto[-1] = new_jobs_auto[-1].replace('JOB_FILE', job_auto_file)
    with open(os.path.join(cali_path, job_auto_file), 'w') as f:
        f.write(new_jobs_auto[-1])
