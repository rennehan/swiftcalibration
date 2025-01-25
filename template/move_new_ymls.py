import os
import numpy as np
from glob import glob
from pathlib import Path

parameter_files = [Path(x) for x in glob("./original_ymls/*.yml")]
parameter_filenames = {filename.stem: filename for filename in parameter_files}

cali_dir = './calibrations'

new_ymls = []
for v,k in enumerate(parameter_filenames):
    cali_name = 'cali_%04d' % int(k)
    cali_path = os.path.join(cali_dir, cali_name)
    yml_file = 'cali_%04d.yml' % int(k)
    os.system('cp ./original_ymls/%d.yml %s' % (int(k), os.path.join(cali_path, yml_file)))
