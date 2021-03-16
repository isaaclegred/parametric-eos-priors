#!/bin/bash

MAX_NUM_SAMPLES=100000

#Sometimes need this
# --copy-column  "logweight_Romani_J1810" 


# These should be the only things that *need* to change
NONPAR_EOS_DIR="/home/philippe.landry/gpr-eos-stacking/EoS/mrgagn"
PARAMETRIC_EOS_DIR="/home/isaac.legred/parametric-eos-priors/good_gaussian_eos_draw_spectral/"
EOS_PER_DIR="1000 100"


# This may need to change if you change the overall convention
TAGS="/corrected_nonparametric_post.csv /corrected_parametric_post.csv"
# This is kinda a silly thing to do, but it will work for now
# Hardcode this if it's easier than worrying about where
# you call functions from
cd ../..
POST_DIR_5=$(pwd)
cd Utils/Plotting



EOS_DIR_TAGS="$NONPAR_EOS_DIR $PARAMETRIC_EOS_DIR"
EOS_COUNT_ARR=($EOS_PER_DIR)
counter=0
for TAG in $TAGS
do
    EOS_DIRS=($EOS_DIR_TAGS)
    EOS_DIR=${EOS_DIRS[$counter]}
    INPATH=$POST_DIR_5$TAG
    echo $INPATH
    OUTPATH="new"$INPATH
    EOS_NUM_PER_DIR=${EOS_COUNT_ARR[$counter]}
    counter=$((counter+1))
    # look up pressures

    # Get the max masses
    $(which weigh-samples) \
        $INPATH \
	$INPATH \
        $OUTPATH \
        Mmax "pressurec2(baryon_density=5.6e+14)" "pressurec2(baryon_density=1.68e+15)" "R(M=1.4)" "Lambda(M=1.4)" \
        --weight-column "logweight_total" \
        --weight-column-is-log "logweight-total" \
        --eos-column eos \
	--Verbose
done
