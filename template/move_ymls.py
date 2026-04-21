import os
# import numpy as np
from glob import glob
from pathlib import Path

params_path = f'/scratch/aspadawe/sims/HyenasC/L1/SimbaC_L1_Calibration/trillium/gizmo-hyenasc-l1-correct_jet-vary_params/cali_params/32_sims-3_params/'

# parameter_files = [Path(x) for x in glob("./original_ymls/*.yml")]
parameter_files = [Path(x) for x in glob(os.path.join(params_path, '*.yml'))]
parameter_filenames = {filename.stem: filename for filename in parameter_files}

# cali_dir = './calibrations'
cali_dir = '/scratch/aspadawe/snapshots/HyenasC/L1/SimbaC_L1_Calibration/halo_3224-correct_jet-calibrations/32_sims-3_params/'

new_ymls = []
for v,k in enumerate(parameter_filenames):
    cali_name = 'cali_%04d' % int(k)
    cali_path = os.path.join(cali_dir, cali_name)
    yml_file = 'cali_%04d.yml' % int(k)
    # os.system('cp ./original_ymls/%d.yml %s' % (int(k), os.path.join(cali_path, yml_file)))
    # os.system('cp %s %s' % (parameter_filenames[k], os.path.join(cali_path, yml_file)))
    os.system(f'cp {parameter_filenames[k]} {os.path.join(cali_path, yml_file)}')
