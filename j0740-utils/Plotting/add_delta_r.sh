#!/bin/bash
TAGS="all_miller all_riley no_j0740"
MAX_NUM_SAMPLES=400000
cd ../..                                                                                       
MAIN_DIR=$(pwd)                                                                                
cd Utils/Plotting

for PRETAG in $TAGS
do
    TAG="/corrected_"$PRETAG"_post.csv"
    samples-calculus \
        $MAIN_DIR"/"$TAG \
        $MAIN_DIR"/"$TAG \
        subtract \
        "R(M=1.4)"\
        "R(M=2.0)" \
        --new-column "DeltaR(2.0-1.4)"\
        --overwrite 
done
