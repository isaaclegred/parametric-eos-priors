#!/bin/bash
counter=$1
echo $counter
MAX_NUM_SAMPLES=350000

#Sometimes need this
# --copy-column  "logweight_Romani_J1810" 


# These should be the only things that *need* to change
NONPAR_EOS_DIR="/home/philippe.landry/nseos/eos/gp/mrgagn/"
NONPAR_EOS_DIR_cs2c2="/home/isaac.legred/local_mrgagn_big_with_cs2c2/"
PARAMETRIC_EOS_DIR="/home/isaac.legred/parametric-eos-priors/good_gaussian_eos_draw_spectral/"
EOS_PER_DIR="1000 1000 1000 1000 1000"

# This may need "/home/isaac.legred/CompParNonParRuns/Plots/All/corrected_nonparametric_post.csv"to change if you change the overall convention
TAGS="all_miller"
echo $TAGS
# This is kinda a silly thing to do, but it will work for now
# Hardcode this if it's easier than worrying about where
# you call functions from
cd ../..
MAIN_DIR=$(pwd)
cd Utils/Plotting

#TAGS="$TAGS qmc-hadagn-0.5rho_nuc-ingo-bps-3e11"
#TAGS="$TAGS qmc-hadagn-0.6rho_nuc-ingo-bps-3e11"
#TAGS="$TAGS qmc-hadagn-0.7rho_nuc-ingo-bps-3e11"
#TAGS="$TAGS qmc-hadagn-0.8rho_nuc-ingo-bps-3e11"
#TAGS="$TAGS qmc-hadagn-0.9rho_nuc-ingo-bps-3e11"
#TAGS="$TAGS qmc-hadagn-1.0rho_nuc-ingo-bps-3e11"
#TAGS="$TAGS qmc-hadagn-1.1rho_nuc-ingo-bps-3e11"

#TAGS="$TAGS mbpt-hadagn-0.5rho_nuc-ingo-bps-3e11"
#TAGS="$TAGS mbpt-hadagn-0.6rho_nuc-ingo-bps-3e11"
#TAGS="$TAGS mbpt-hadagn-0.7rho_nuc-ingo-bps-3e11"
#TAGS="$TAGS mbpt-hadagn-0.8rho_nuc-ingo-bps-3e11"
#TAGS="$TAGS mbpt-hadagn-0.9rho_nuc-ingo-bps-3e11"
#TAGS="$TAGS mbpt-hadagn-1.0rho_nuc-ingo-bps-3e11"
#TAGS="$TAGS mbpt-hadagn-1.1rho_nuc-ingo-bps-3e11"
#TAGS="$TAGS mbpt-hadagn-1.2rho_nuc-ingo-bps-3e11"
#TAGS="$TAGS mbpt-hadagn-1.3rho_nuc-ingo-bps-3e11"
#TAGS="$TAGS mbpt-hadagn-1.4rho_nuc-ingo-bps-3e11"

#TAGS="$TAGS mrgagn-ingo-bps-1e10"



EOS_DIR_TAGS="$NONPAR_EOS_DIR $NONPAR_EOS_DIR $NONPAR_EOS_DIR $NONPAR_EOS_DIR $NONPAR_EOS_DIR"
EOS_CS2C2_DIR_TAGS="$NONPAR_EOS_DIR_cs2c2 $NONPAR_EOS_DIR_cs2c2 $NONPAR_EOS_DIR_cs2c2 $NONPAR_EOS_DIR_cs2c2 $NONPAR_EOS_DIR_cs2c2"
EOS_COUNT_ARR=($EOS_PER_DIR)

PRETAGS=($TAGS)
PRETAG=${PRETAGS[$counter]}

TAG="/corrected_"$PRETAG"_post.csv"
OUTTAG="/corrected_"$PRETAG"_post_test.csv"
EOS_DIRS=($EOS_DIR_TAGS)
EOS_DIR=${EOS_DIRS[$counter]}
EOS_DIRS_CS2=($EOS_CS2C2_DIR_TAGS)
EOS_DIR_cs2c2=${EOS_DIRS_CS2[$counter]}
INPATH=$MAIN_DIR$TAG
OUTPATH=$MAIN_DIR$OUTTAG
EOS_NUM_PER_DIR=${EOS_COUNT_ARR[$counter]}


$(which process2count) \
    $INPATH\
    $OUTPATH \
    branch\
    --copy-column  "logweight_total" \
    --copy-column "Mmax" \
    --copy-column "pressurec2(baryon_density=2.8e+14)" \
    --copy-column "pressurec2(baryon_density=5.6e+14)"\
    --copy-column "pressurec2(baryon_density=1.68e+15)"\
    --copy-column "pressurec2(baryon_density=1.96e+15)"\
    --copy-column "energy_densityc2(baryon_density=2.8e+14)"\
    --copy-column "energy_densityc2(baryon_density=5.6e+14)"\
    --copy-column "energy_densityc2(baryon_density=1.68e+15)" \
    --copy-column "energy_densityc2(baryon_density=1.96e+15)" \
    --copy-column "R(M=1.4)" \
    --copy-column "R(M=2.0)" \
    --copy-column "Lambda(M=1.4)" \
    --copy-column "Lambda(M=2.0)" \
    --copy-column "DeltaR(2.0-1.4)" \
    --copy-column "rhoc(M@Mmax)" \
    --copy-column "cs2c2max"\
    --copy-column "baryon_density(cs2c2@cs2c2max)"\
    --copy-column "pressurec2(cs2c2@cs2c2max)"\
    --eos-column eos \
    --eos-dir ${EOS_DIR_cs2c2} \
    --eos-num-per-dir ${EOS_NUM_PER_DIR} \
    --eos-basename 'macro-draw-%(draw)06d-branches.csv' \
    --max-num-samples $MAX_NUM_SAMPLES \
    --Verbose

