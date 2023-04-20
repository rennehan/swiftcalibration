import os
import numpy as np
from glob import glob
from pathlib import Path

data_dir = '/home/b/babul/rennehan/work/swift'
files_to_link = ['swift', 'yieldtables', 'photometry', 'output_list_cali.txt', 'UV_dust1_CR1_G1_shield1.hdf5', 'Jenny_N128L12.dat.hdf5', 'velo.cfg']

with open('job.sh', 'r') as f:
    job_data = f.read()

parameter_files = [Path(x) for x in glob("./original_ymls/*.yml")]
parameter_filenames = {filename.stem: filename for filename in parameter_files}

new_jobs = []
new_ymls = []
for v,k in enumerate(parameter_filenames):
    cali_name = 'cali_%04d' % int(k)
    yml_file = 'swimba_%04d.yml' % int(k)
    job_file = 'job_%04d.sh' % int(k)
    os.makedirs(cali_name, mode = 0o755, exist_ok = True)
    for file_to_link in files_to_link:
        os.system('ln -s %s/%s %s/%s' % (data_dir, file_to_link, cali_name, file_to_link))

    new_jobs.append(job_data.replace('JOB_NAME', cali_name))
    new_jobs[-1] = new_jobs[-1].replace('YML_FILE', yml_file)
    with open('%s/%s' % (cali_name, job_file), 'w') as f:
        f.write(new_jobs[-1])

    os.system('cp ./original_ymls/%d.yml %s/%s' % (int(k), cali_name, yml_file))

