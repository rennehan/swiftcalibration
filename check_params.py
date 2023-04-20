import os
import numpy as np
from glob import glob
from pathlib import Path
import yaml

parameter_files = [Path(x) for x in glob("./original_ymls/*.yml")]
parameter_filenames = {filename.stem: filename for filename in parameter_files}

keys = ["torque_accretion_norm",
        "slim_disk_wind_speed_km_s",
        "adaf_f_accretion",
        "quasar_wind_speed_km_s",
        "jet_velocity",
        "adaf_coupling"]

values = [[], [], [], [], [], []]

for i,k in enumerate(parameter_filenames):
    data_dict = yaml.load(open('./original_ymls/%d.yml' % int(k), 'r'),
                          Loader=yaml.FullLoader)
    for j,key in enumerate(keys):
        if key == 'adaf_f_accretion':
            print(k, data_dict['YAMAGN'][key], data_dict['YAMAGN']['jet_velocity'])
        values[j].append(data_dict['YAMAGN'][key])

