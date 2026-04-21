import os
import numpy as np
from glob import glob
from pathlib import Path
import copy



print("\nGenerating job files...\n")


## Specify paths and parameters for generating calibrations; modify as needed for your calibration parameters and model
job_template = '/path/to/job_template.sh'  # template for job script to run each calibration; should contain placeholders for job name, number of jobs, number of nodes, and parameter file name
job_restart_template = '/path/to/job_restart_template.sh'  # template for job script to restart each calibration; should contain placeholders for job name, number of jobs, number of nodes, and parameter file name

cali_dir = '/path/to/calibration/parameter/files/'  # same as output_dir in design_calibrations.py
output_dir = '/path/to/output/calibration/directory/'  # same as output_dir in generate_calibrations.py
param_type = 'tex'

num_nodes = 3  # number of nodes to request for each slurm job
num_subjobs = 15  # number of sub-jobs in each slurm job array
job_name_specifier = '32_sims-3_params'  # string to include in job name to specify calibration setup; modify as needed for your calibration parameters and model



## Check if output directory exists; if not, create it
## Note: this should already be done in generate_calibrations.py, so this is just a safety check
if not os.path.isdir(output_dir):
    print(f"\tOutput directory {output_dir} does not exist; creating it.")
    os.makedirs(output_dir, mode = 0o755, exist_ok = True)



# parameter_files = [Path(x) for x in glob("./original_ymls/*.yml")]
# parameter_files = [Path(x) for x in glob(f"{cali_dir}/*.{param_type}")]
# parameter_filenames = {filename.stem: filename for filename in parameter_files}

# ## Find parameter files in cali_dir with specified type (yml or tex) and generate calibration setups for each parameter file
# parameter_files = [Path(x) for x in glob(os.path.join(cali_dir, f"*.{param_type}"))]  # look for parameter files in cali_dir with format *.yml or *.tex
# parameter_filenames = {filename.stem:filename for filename in parameter_files}
# print(f"\tFound {len(parameter_files)} parameter files in {cali_dir} with type {param_type}")


## Set up dictionary with calibration info
## Change param_file, job_file, and job_restart_file as needed for your calibration setup
cali_subdirs = [Path(x) for x in glob(os.path.join(output_dir, 'cali_*')) if os.path.isdir(x)]
cali_dict = {
    cali_dir.stem:{
        'cali_dir':cali_dir,
        'param_file':f'params.{param_type}',
        'job_file':f'job.sh',
        'job_restart_file':f'job_restart.sh',
    } 
    for cali_dir in cali_subdirs
}



## Read in job templates
with open(job_template, 'r') as f:
    job_data = f.read()

with open(job_restart_template, 'r') as f:
    job_restart_data = f.read()


## Generate job files for each calibration setup
# new_jobs = []
# new_jobs_restart = []
for i, (cali_name, cali_info) in enumerate(cali_dict.items()):
    print(f"\tGenerating job files for calibration setup {i+1}/{len(cali_dict)}: {cali_name}...")

    cali_dir = cali_info['cali_dir']
    param_file = cali_info['param_file']
    job_file = cali_info['job_file']
    job_restart_file = cali_info['job_restart_file']


    new_job = copy.deepcopy(job_data)
    new_job = new_job.replace('JOB_NAME', f'{cali_name}-{job_name_specifier}')
    new_job = new_job.replace('NUM_JOBS', str(num_subjobs))
    new_job = new_job.replace('NUMBER_NODES', str(num_nodes))
    new_job = new_job.replace('PARAMETER_FILE', param_file)
    with open(os.path.join(cali_dir, job_file), 'w') as f:
        f.write(new_job)

    new_job_restart = copy.deepcopy(job_restart_data)
    new_job_restart = new_job_restart.replace('JOB_NAME', f'{cali_name}-{job_name_specifier}-restart')
    new_job_restart = new_job_restart.replace('NUM_JOBS', str(num_subjobs))
    new_job_restart = new_job_restart.replace('NUMBER_NODES', str(num_nodes))
    new_job_restart = new_job_restart.replace('PARAMETER_FILE', param_file)
    with open(os.path.join(cali_dir, job_restart_file), 'w') as f:
        f.write(new_job_restart)


    # new_jobs.append(job_data.replace('JOB_NAME', f'{cali_name}-{job_name_specifier}'))
    # new_jobs[-1] = new_jobs[-1].replace('NUM_JOBS', str(num_subjobs))
    # new_jobs[-1] = new_jobs[-1].replace('NUMBER_NODES', str(num_nodes))
    # new_jobs[-1] = new_jobs[-1].replace('PARAMETER_FILE', param_file)
    # with open(os.path.join(cali_dir, job_file), 'w') as f:
    #     f.write(new_jobs[-1])

    # new_jobs_restart.append(job_restart_data.replace('JOB_NAME', f'{cali_name}-{job_name_specifier}-restart'))
    # new_jobs_restart[-1] = new_jobs_restart[-1].replace('NUM_JOBS', str(num_subjobs))
    # new_jobs_restart[-1] = new_jobs_restart[-1].replace('NUMBER_NODES', str(num_nodes))
    # new_jobs_restart[-1] = new_jobs_restart[-1].replace('PARAMETER_FILE', param_file)
    # with open(os.path.join(cali_dir, job_restart_file), 'w') as f:
    #     f.write(new_jobs_restart[-1])


print("\nFinished generating job files.\n")