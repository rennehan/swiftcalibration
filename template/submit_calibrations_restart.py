import os
import numpy as np

num_simulations = 32
calis_to_submit = np.arange(0, num_simulations)

## Can also specify particular calibrations to submit by changing calis_to_submit to a list of calibration numbers; e.g.
# calis_to_submit = [0, 1, 2, 3]
## Or can specify certain calibrations NOT to resubmit, e.g.
#calis_to_submit = np.delete(calis_to_submit, [0,1,2,3])

print(calis_to_submit)

output_dir = '/path/to/output/calibration/directory/'  # same as output_dir in generate_calibrations.py

for cali_to_submit in calis_to_submit:
    ## Remember to change job_restart.sh if used a different job name in generate_jobs.py; e.g. job_restart_0000.sh
    os.system(f'cd {output_dir}/cali_{cali_to_submit:04d} && sbatch job_restart.sh')
