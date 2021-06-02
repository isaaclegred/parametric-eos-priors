#!/bin/bash

### test extracting quantiles
############################
# THIS WILL HAVE TO CHANGE TO BE THE
# ONE COPIED OVER WITH CS2C2
############################\/\/\/\/\/\/\/\/
export NONPAR_EOS_DIR="/home/isaac.legred/local_mrgagn_big_with_cs2c2/"
############################/\/\/\/\/\/\/\/\
EOS_PER_DIR="1000 1000 1000 1000"
TAGS="psr psrgw all_miller no_j0740"

cd ../..
MAIN_DIR=$(pwd)
cd Utils/Plotting


# #echo \
EOS_DIR_TAGS="$NONPAR_EOS_DIR $NONPAR_EOS_DIR $NONPAR_EOS_DIR $NONPAR_EOS_DIR"
EOS_COUNT_ARR=($EOS_PER_DIR)
counter=0

for PRETAG in $TAGS
do
    TAG="/corrected_"$PRETAG"_post.csv"
    EOS_DIRS=($EOS_DIR_TAGS)
    EOS_DIR=${EOS_DIRS[$counter]}
    INPATH=$POST_DIR_5$TAG
    echo $INPATH
    OUTPATH=$INPATH
    EOS_NUM_PER_DIR=${EOS_COUNT_ARR[$counter]}
    counter=$((counter+1))
    
    process2quantiles \
        $MAIN_DIR"/"$TAG \
        $MAIN_DIR"/"$PRETAG"_quantiles_cs2.csv"\
        baryon_density \
        cs2c2 \
        2.8e13 2.8e15 \
        --logcolumn baryon_density \
        --max-num-samples 200000 \
        --weight-column logweight_total \
        --weight-column-is-log logweight_total \
        --eos-column eos \
        --eos-dir $EOS_DIR \
        --eos-num-per-dir $EOS_NUM_PER_DIR \
        --eos-basename 'eos-draw-%(draw)06d.csv' \
        --num-points 100 \
        --Verbose 
done

# Prior is special case
process2quantiles \
    $MAIN_DIR"/corrected_all_miller_post.csv" \
    $MAIN_DIR"/prior_quantiles_cs2.csv"\
    baryon_density \
    cs2c2 \
    2.8e13 2.8e15 \
    --logcolumn baryon_density \
    --max-num-samples 200000 \
    --eos-column eos \
    --eos-dir $NONPAR_EOS_DIR \
    --eos-num-per-dir 1000 \
    --eos-basename 'eos-draw-%(draw)06d.csv' \
    --num-points 100 \
    --Verbose
