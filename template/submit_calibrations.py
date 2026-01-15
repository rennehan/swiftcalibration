import os
import numpy as np

calis_to_submit = np.arange(0, 32)
print(calis_to_submit)

# cali_dir = './calibrations'
output_dir = '/scratch/aspadawe/snapshots/HyenasC/L1/SimbaC_L1_Calibration/halo_3224-correct_jet-calibrations/32_sims-6_params/'

for cali_to_submit in calis_to_submit:
    # os.system('cd %s/cali_%04d && sbatch job_%04d.sh' % (output_dir, cali_to_submit, cali_to_submit))
    os.system(f'cd {output_dir}/cali_{cali_to_submit:04d} && sbatch job_{cali_to_submit:04d}.sh')