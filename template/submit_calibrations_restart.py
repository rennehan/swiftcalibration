import os
import numpy as np

calis_to_submit = np.arange(0, 32)

cali_dir = './calibrations'

for cali_to_submit in calis_to_submit:
    os.system('cd %s/cali_%04d && sbatch job_restart_%04d.sh' % (cali_dir, cali_to_submit, cali_to_submit))
