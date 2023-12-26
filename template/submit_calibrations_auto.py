import os
import numpy as np

calis_to_submit = np.arange(0, 32)

cali_dir = './calibrations'

for cali_to_submit in calis_to_submit:
    ## 0 for starting with no restart files initially, 1 for starting from restart files initially
    os.system('cd %s/cali_%04d && sbatch job_auto_%04d.sh 1' % (cali_dir, cali_to_submit, cali_to_submit))
