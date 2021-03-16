# Create a copy of a posterior file but with an additional column
# which is some prefined function applied to already exisiting columns 
import numpy as np
import argparse

# example function (should be vectorized otherwise will not work)
def compactness_lam(lam):
    # Universal relation, eq. 78 in YY https://arxiv.org/pdf/1608.02582.pdf
    a_vals= [.360, -.0355, .000705]
    output =  np.sum([a_vals[k]*lam**k for k in range(3)], 0)
    return output
parser = argparse.ArgumentParser(description="Get File Name")

parser.add_argument('--posterior-file-path', dest="source_path", default="./total_post.csv",
                    help = "The path to post.csv file that contains columns that can be copied \
                    into the prior file ")

parser.add_argument('--output', dest="output", default="./total_prior.csv",
                    help = "Where to create the output (prior) csv file ")
parser.add_argument('--new-cov-name', dest="new_cov_name", default="new_cov",
                    help = "The name of the new covariate in the file")
parser.add_argument('--add-parameter', dest="parameters", default=[], action='append', help = "add a parameter to evaluate the function")

args = parser.parse_args()

source_path = args.source_path
output_path = args.output

to_call = compactness_lam
parameters = args.parameters
f = open(source_path)
source_header = f.readline()
header_list= source_header.split(",")
new_var_name = args.new_cov_name
f.close()

param_indices = []
for i, param in enumerate(parameters):
    param_indices.append(header_list.index(param))
print(len(param_indices))

source_data = np.loadtxt(source_path,skiprows=1,delimiter=",")        
num_points = len(source_data[:,0])
# Make a table of all the variable values we could possibly need
var_vals = []
for i, param_index in enumerate(param_indices):
    var_vals.append(source_data[:, param_index])

new_vals = to_call(*var_vals)
new_vals = np.reshape(new_vals, (len(new_vals),1))
new_header = source_header[:-1] + "," + new_var_name + "\n"
new_data=np.concatenate([source_data, new_vals], 1)
np.savetxt(output_path, source_data,  header = new_header, comments="", delimiter=",")
