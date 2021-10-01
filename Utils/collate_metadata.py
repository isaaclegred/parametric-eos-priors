import numpy as np
import pandas as pd

num_dirs = 2000
eos_per_dir = 100

total_eos = num_dirs * eos_per_dir

eos_dir = "/home/isaac.legred/parametric-eos-priors/eos_draws/production_eos_draw_sos"
num_params = 6

# We need one column for each param, +1  for the EoS column
collated_data = np.ndarray((total_eos, num_params + 1))
header=""
for dir_index in range(num_dirs):
    file_name = eos_dir + "/DRAWmod" + str(eos_per_dir) + "-%06d" %(dir_index)+"/eos_metadata-" + "%06d" %(dir_index) + ".csv"
    if dir_index == 0:
        file1 = open(file_name)
        header = file1.readline()[:-1]
    local_metadata = np.loadtxt(file_name, skiprows=1, comments="", delimiter=",")
    collated_data[dir_index*eos_per_dir : (dir_index + 1)*eos_per_dir] = local_metadata

np.savetxt(eos_dir + "/all_eos_metadata.csv", collated_data, header=header, comments="", delimiter=",")
