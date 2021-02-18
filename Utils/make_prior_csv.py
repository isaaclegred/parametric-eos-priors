# Extract the relevant columns from posterior csv to create an analogous prior csv that can be
# used in the plotting functions.  
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Get File Name")
parser.add_argument('--posterior-file-path', dest="source_path", default="./total_post.csv",
                    help = "The path to post.csv file that contains columns that can be copied \
                    into the prior file ")

parser.add_argument('--output', dest="output", default="./total_prior.csv",
                    help = "Where to create the output (prior) csv file ")
parser.add_argument('--num-weight-columns', dest="num_weights", default=-1, type=int,
                    help="Number of logweight columns to copy")
args = parser.parse_args()

source_path = args.source_path
output_path = args.output
num_weights = args.num_weights
f = open(source_path)
source_header = f.readline()
f.close() 
source_data = np.loadtxt(source_path,skiprows=1,delimiter=",")
eos_indices = source_data[:,0]
if (num_weights == -1):
    prior_data = np.ones_like(source_data)
    prior_data[:,0] = eos_indices
else:
    weight_reps  = np.ones_like(source_data[:,1:num_weights+1])
    prior_data = source_data
    prior_data[:,0] = eos_indices
    prior_data[:,1:num_weights+1] = weight_reps
    
np.savetxt(output_path, prior_data,  header = source_header, comments="", delimiter=",")
