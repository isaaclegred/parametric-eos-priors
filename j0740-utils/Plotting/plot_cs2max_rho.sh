#!/bin/bash

### make plots of nuclear parameters alongside EoS params and macroscopic observables

# E_beta,E_PNM,x,S0,L,Ksym,eps_total,eps_electron,mu_electron,EoS,Mmax,numbranches,logweight_NICER,logweight_GW170817,logweight_GW190425,logweight_J0348,logweight_J0740,logweight_J1614,pressurec2(baryon_density=1.4e+14),pressurec2(baryon_density=2.8e+14),pressurec2(baryon_density=4.2e+14),pressurec2(baryon_density=5.6e+14),pressurec2(baryon_density=8.4e+14),pressurec2(baryon_density=1.12e+15),energy_densityc2(baryon_density=1.4e+14),energy_densityc2(baryon_density=2.8e+14),energy_densityc2(baryon_density=4.2e+14),energy_densityc2(baryon_density=5.6e+14),energy_densityc2(baryon_density=8.4e+14),energy_densityc2(baryon_density=1.12e+15),R(M=1.4),R(M=1.6),Lambda(M=1.4),Lambda(M=1.6)

#-------------------------------------------------

MAX_NUM_SAMPLES=100000
NUM_POINTS=101
#---

# plot different posteriors for the same set of prior draws
TAGS="/home/philippe.landry/nseos/eos/gp/mrgagn/"

for TAG in $TAGS
do

    echo "----------"
    echo "processing $TAG"
    echo "----------"
    # Again this is silly but will make things easier
    # This script should only be called by condor, and in particular with
    # the submit file that requires you to be in the correct directory anways
    # this shouldn't cause big problems for now
    cd ../..
    MAIN_DIR=$(pwd)
    cd Utils/Plotting
    
    ALL_SAMPS=$MAIN_DIR"/corrected_all_miller_post.csv"
    PRIOR_SAMPS=$MAIN_DIR"/corrected_no_j0740_post.csv"
    NO_NICER_SAMPS=$MAIN_DIR"/corrected_no_j0740_post.csv"

    PRIOR="PRIOR"
    PSR="PSR"
    PSRGW="PSR+GW"
    ALL="PSR+GW+J0030+J0740"
    NO_NICER="PSR+GW+J0030"
    #ALL="Marginalized-Composition"
    SINGLE="1-Branch"
    MULTI="2+-Branches"



    kde-corner-samples \
        'cs2c2max' 'baryon_density(cs2c2@cs2c2max)'  \
        --column-label 'cs2c2max' '$c^2_{s,\mathrm{max}}/c^2$' \
        --column-label 'baryon_density(cs2c2@cs2c2max)' '$\log_{10}\left(\rho(c^2_{s,\mathrm{max}}/c^2)\ [\mathrm{g}/\mathrm{cm}]\right)$' \
        -s $PRIOR $PRIOR_SAMPS \
        -s $NO_NICER $NO_NICER_SAMPS \
        -s $ALL $ALL_SAMPS \
        --column-range 'baryon_density(cs2c2@cs2c2max)' 14.5 15.5\
        --logcolumn 'baryon_density(cs2c2@cs2c2max)'\
        --column-multiplier 'baryon_density(cs2c2@cs2c2max)' .4343\
        --weight-column $ALL logweight_total\
        --weight-column-is-log $ALL logweight_total\
        --weight-column $NO_NICER logweight_total\
        --weight-column-is-log $NO_NICER logweight_total\
        --color $ALL "b" \
        --color $PRIOR "k" \
        --color $NO_NICER "c" \
        -B 'cs2c2max' .330 .335\
        -B 'baryon_density(cs2c2@cs2c2max)' 5.57e14 5.63e14\
        -B 'baryon_density(cs2c2@cs2c2max)' 1.65e15 1.71e15\
        --output-dir $MAIN_DIR \
        --tag "covariates-cs2c2max-rhoc_at_cs2max" \
        --num-proc 4 \
        --no-scatter\
        --num-points $NUM_POINTS\
        --level 0.9 \
        --legend \
        --grid \
        --rotate-xticklabels 90 \
        --figwidth 5.5 \
        --figheight 5.5 \
        --Verbose         
done
