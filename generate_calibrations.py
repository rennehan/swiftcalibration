import os
import numpy as np
from glob import glob
from pathlib import Path
import argparse as ap

# System = niagara, rusty
# data_dir = Folder with all of the data for
#            running the simulation. Need to
#            link these files
# ics_file = The HDF5 file with the ICs.
parser = ap.ArgumentParser()
parser.add_argument("system")
parser.add_argument("data_dir")
parser.add_argument("ics_file")
args = parser.parse_args()

files_to_link = ["swift", 
                 "yieldtables", 
                 "photometry", 
                 "output_list_cali.txt", 
                 "UV_dust1_CR1_G1_shield1.hdf5", 
                 args.ics_file]

with open("job_%s.sh" % args.system, "r") as f:
    job_data = f.read()

parameter_files = [Path(x) for x in glob("./original_ymls/*.yml")]
parameter_filenames = {filename.stem: filename for filename in parameter_files}

new_jobs = []
new_ymls = []
for v,k in enumerate(parameter_filenames):
    cali_name = "cali_%04d" % int(k)
    yml_file = "swimba_%04d.yml" % int(k)
    job_file = "job_%04d.sh" % int(k)
    os.makedirs(cali_name, mode = 0o755, exist_ok = True)
    for file_to_link in files_to_link:
        os.system("ln -s %s/%s %s/%s" % (args.data_dir, file_to_link, cali_name, file_to_link))

    new_jobs.append(job_data.replace("JOB_NAME", cali_name))
    new_jobs[-1] = new_jobs[-1].replace("YML_FILE", yml_file)
    with open("%s/%s" % (cali_name, job_file), "w") as f:
        f.write(new_jobs[-1])

    os.system("cp ./original_ymls/%d.yml %s/%s" % (int(k), cali_name, yml_file))

