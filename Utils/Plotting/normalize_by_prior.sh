#!/bin/bash
MAX_NUM_SAMPLES=100000

#Sometimes need this
# --copy-column  "logweight_Romani_J1810" 


# This may need to change if you change the overall convention
TAGS="/corrected_nonparametric_post.csv /corrected_parametric_post.csv"
BANDWIDTH_TAGS=".09 .07"
COLUMN_NAMES="Mmax"
# This is kinda a silly thing to do, but it will work for now
# Hardcode this if it's easier than worrying about where
# you call functions from
cd ../..
POST_DIR_5=$(pwd)
cd Utils/Plotting




counter=0
for TAG in $TAGS
do
    BANDWIDTHS=($BANDWIDTH_TAGS)
    BANDWIDTH=${BANDWIDTHS[$counter]}
    INPATH=$POST_DIR_5$TAG
    echo $INPATH
    OUTPATH=$POST_DIR_5"normalized_by_prior"$TAG
    EOS_NUM_PER_DIR=${EOS_COUNT_ARR[$counter]}
    counter=$((counter+1))
    # look up pressures

    # Get the max masses
    $(which weigh-samples) \
        $INPATH \
	$INPATH \
        $INPATH \
        $COLUMN_NAMES \
        --bandwidth $BANDWIDTH \
        --output-weight-column $COLUMN_NAMES"_prior_kde"\
        -v
done
