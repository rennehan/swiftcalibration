import os
import numpy as np
from glob import glob
from pathlib import Path


job_template = '/scratch/aspadawe/sims/HyenasC/L1/SimbaC_L1_Calibration/trillium/gizmo-hyenasc-l1-correct_jet-vary_params/job_array-template.sh'
job_restart_template = '/scratch/aspadawe/sims/HyenasC/L1/SimbaC_L1_Calibration/trillium/gizmo-hyenasc-l1-correct_jet-vary_params/job_array_restart-template.sh'

cali_dir = '/scratch/aspadawe/sims/HyenasC/L1/SimbaC_L1_Calibration/trillium/gizmo-hyenasc-l1-correct_jet-vary_params/cali_params/128_sims-6_params/'
param_type = 'tex'
output_dir = '/scratch/aspadawe/snapshots/HyenasC/L1/SimbaC_L1_Calibration/halo_2205-correct_jet-calibrations/128_sims-6_params/'

num_nodes = 6
num_jobs = 5


# with open('job.sh', 'r') as f:
#     job_data = f.read()

# with open('job_restart.sh', 'r') as f:
#     job_restart_data = f.read()

# with open('job_mpi.sh', 'r') as f:
#     job_mpi_data = f.read()

# with open('job_restart_mpi.sh', 'r') as f:
#     job_restart_mpi_data = f.read()

# with open('job_auto.sh', 'r') as f:
#     job_auto_data = f.read()


# parameter_files = [Path(x) for x in glob("./original_ymls/*.yml")]
parameter_files = [Path(x) for x in glob(f"{cali_dir}/*.{param_type}")]
parameter_filenames = {filename.stem: filename for filename in parameter_files}

with open(job_template, 'r') as f:
    job_data = f.read()

with open(job_restart_template, 'r') as f:
    job_restart_data = f.read()

# cali_dir = './calibrations'
if not os.path.isdir(output_dir):
    os.makedirs(output_dir, mode = 0o755, exist_ok = True)

new_jobs = []
new_jobs_restart = []
# new_jobs_mpi = []
# new_jobs_restart_mpi = []
# new_jobs_auto = []
# new_ymls = []
for v,k in enumerate(parameter_filenames):
    # cali_name = 'cali_%04d' % int(k)
    # cali_path = os.path.join(output_dir, cali_name)
    # param_file = 'cali_%04d.yml' % int(k)
    cali_name = f'cali_{int(k):04d}'
    cali_path = os.path.join(output_dir, cali_name)
    param_file = f'cali_{int(k):04d}.{param_type}'

    # job_file = 'job_%04d.sh' % int(k)
    job_file = f'job_{int(k):04d}.sh'
    job_restart_file = f'job_restart_{int(k):04d}.sh'
    # job_mpi_file = 'job_mpi_%04d.sh' % int(k)
    # job_restart_mpi_file = 'job_restart_mpi_%04d.sh' % int(k)
    # job_auto_file = 'job_auto_%04d.sh' % int(k)


    new_jobs.append(job_data.replace('JOB_NAME', cali_name))
    new_jobs[-1] = new_jobs[-1].replace('NUM_JOBS', str(num_jobs))
    new_jobs[-1] = new_jobs[-1].replace('NUMBER_NODES', str(num_nodes))
    new_jobs[-1] = new_jobs[-1].replace('PARAMETER_FILE', param_file)
    with open(os.path.join(cali_path, job_file), 'w') as f:
        f.write(new_jobs[-1])

    new_jobs_restart.append(job_restart_data.replace('JOB_NAME', cali_name))
    new_jobs_restart[-1] = new_jobs_restart[-1].replace('NUM_JOBS', str(num_jobs))
    new_jobs_restart[-1] = new_jobs_restart[-1].replace('NUMBER_NODES', str(num_nodes))
    new_jobs_restart[-1] = new_jobs_restart[-1].replace('PARAMETER_FILE', param_file)
    with open(os.path.join(cali_path, job_restart_file), 'w') as f:
        f.write(new_jobs_restart[-1])

    # new_jobs_mpi.append(job_mpi_data.replace('JOB_NAME', cali_name))
    # new_jobs_mpi[-1] = new_jobs_mpi[-1].replace('YML_FILE', yml_file)
    # new_jobs_mpi[-1] = new_jobs_mpi[-1].replace('JOB_FILE', job_mpi_file)
    # with open(os.path.join(cali_path, job_mpi_file), 'w') as f:
    #     f.write(new_jobs_mpi[-1])

    # new_jobs_restart_mpi.append(job_restart_mpi_data.replace('JOB_NAME', cali_name))
    # new_jobs_restart_mpi[-1] = new_jobs_restart_mpi[-1].replace('YML_FILE', yml_file)
    # new_jobs_restart_mpi[-1] = new_jobs_restart_mpi[-1].replace('JOB_FILE', job_restart_mpi_file)
    # with open(os.path.join(cali_path, job_restart_mpi_file), 'w') as f:
    #     f.write(new_jobs_restart_mpi[-1])
    
    # new_jobs_auto.append(job_auto_data.replace('JOB_NAME', cali_name))
    # new_jobs_auto[-1] = new_jobs_auto[-1].replace('YML_FILE', yml_file)
    # new_jobs_auto[-1] = new_jobs_auto[-1].replace('JOB_FILE', job_auto_file)
    # with open(os.path.join(cali_path, job_auto_file), 'w') as f:
    #     f.write(new_jobs_auto[-1])
