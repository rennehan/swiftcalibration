import os
import numpy as np

num_simulations = 32
calis_to_submit = np.arange(0, num_simulations)
print(calis_to_submit)

output_dir = '/path/to/output/calibration/directory/'  # same as output_dir in generate_calibrations.py
# job_name = 'job.sh'

for cali_to_submit in calis_to_submit:
    ## Remember to change job.sh if used a different job name in generate_jobs.py; e.g. job_0000.sh
    os.system(f'cd {output_dir}/cali_{cali_to_submit:04d} && sbatch job.sh')