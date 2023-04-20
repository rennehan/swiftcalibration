import os
import numpy as np

calis_to_submit = np.arange(0, 64)

for cali_to_submit in calis_to_submit:
    os.system('cd cali_%04d && sbatch job_%04d.sh' % (cali_to_submit, cali_to_submit))

