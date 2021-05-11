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
header_list= source_header.split(",")
f.close()
new_header=header_list[0] + ","

for index in range(1, num_weights+1):
    current_event = header_list[index]
    if "total" in current_event:
        new_header += "logweight_total,"
    if "Antoniadis_J0348" in current_event:
        new_header += "logweight_Antoniadis_J0348,"
    if  "LVC_GW170817_PhenomPNRTlo" in current_event:
        new_header += "logweight_LVC_GW170817_PhenomPNRTlo,"
    if  "Miller_J0030_threespot" in current_event:
        new_header += "logweight_Miller_J0030_threespot,"
    if "LVC_GW190425_PhenomPNRThi"in current_event:
        new_header += "logweight_LVC_GW190425_PhenomPNRThi,"
    if "Cromartie_J0740" in current_event:
        new_header += "logweight_Cromartie_J0740,"
    if "Romani_J1810" in current_event:
        new_header += "logweight_Romani_J1810,"
    if "Amsterdam_J0740" in current_event:
        new_header += "logweight_Amsterdam_J0740,"        
if new_header[-1] == ",":
    new_header = new_header[:-1]

        
source_data = np.loadtxt(source_path,skiprows=1,delimiter=",")    
np.savetxt(output_path, source_data,  header = new_header, comments="", delimiter=",")
