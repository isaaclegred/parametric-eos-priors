# Use data from an exisitng CSV to apply a function on certain of columns 
# And then save the used columns, and computed columns to a new csv file
import numpy as np
import argparse

# example function (should be vectorized otherwise will not work)
def R_lam(lam):
    # Universal relation, eq. 78 in YY https://arxiv.org/pdf/1608.02582.pdf
    a_vals= [.360, -.0355, .000705]
    C =  np.sum([a_vals[k]*np.log(lam)**k for k in range(3)], 0)
    # Definition of compacnes, sun is 1.5 km in geometrized units, and we want
    # 1.4 solar mass benchmark
    Rlam = 1.4*1.5/C 
    return Rlam
# example function (should be vectorized otherwise will not work)
def R_lam(lam):
    # Universal relation, eq. (9) https://arxiv.org/pdf/1903.03909.pdf
    a_vals= [.3617, -.0348, .0006194]
    C =  np.sum([a_vals[k]*np.log(lam)**k for k in range(3)], 0)
    # Definition of compacnes, sun is 1.5 km in geometrized units, and we want
    # 1.4 solar mass benchmark
    Rlam = 1.4*1.5/C
    return Rlam


parser = argparse.ArgumentParser(description="Get File Name")

parser.add_argument('--posterior-file-path', dest="source_path", default="./total_post.csv",
                    help = "The path to post.csv whose columns will be used to compute a thing ")

parser.add_argument('--output', dest="output", default="./total_prior.csv",
                    help = "Where to create the reference  csv file ")
parser.add_argument('--new-cov-name', dest="new_cov_name", default="new_cov",
                    help = "The name of the new covariate in the file")
parser.add_argument('--add-parameter', dest="parameters", default=[], action='append', help = "add a parameter to evaluate the function")

args = parser.parse_args()

source_path = args.source_path
output_path = args.output
to_eval_on = np.linspace(0,2200,100)
def prepare_header_from_list(l):
    s = ""
    for elt in l:
        s += elt + ","
    return s
to_call = R_lam
parameters = args.parameters
f = open(source_path)
source_header = f.readline()
header_list= source_header.split(",")
new_var_name = args.new_cov_name
f.close()

param_indices = []
for i, param in enumerate(parameters):
    param_indices.append(header_list.index(param))


source_data = np.loadtxt(source_path,skiprows=1,delimiter=",")        
num_points = len(source_data[:,0])
# Make a table of all the variable values we could possibly need
var_vals = []
for i, param_index in enumerate(param_indices):
    var_vals.append(source_data[:, param_index])

new_vals = to_call(to_eval_on)
new_vals = np.reshape(new_vals, (len(new_vals),1))
print(new_vals.shape)
varr_vals = np.reshape(to_eval_on, (len(new_vals),1))
print(varr_vals.shape)

new_header = prepare_header_from_list(parameters) + new_var_name 
new_data=np.transpose(np.concatenate([np.transpose(varr_vals), np.transpose(new_vals)]))
print(new_data.shape)
np.savetxt(output_path, new_data,  header = new_header, comments="", delimiter=",")
