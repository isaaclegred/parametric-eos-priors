#!/bin/bash

TAGS="parametric nonparametric"
COLUMNS="Mmax"
RANGE=".01 .8"
WRITETAG=$COLUMNS"_prior"

for TAG in $TAGS
do
    INPATH="corrected_"$TAG"_post.csv"
    COLUMNS="Mmax"
    echo $COLUMNS
    optimize-bandwidth \
        $INPATH \
        $COLUMNS \
        $RANGE \
        --num-withheld 100 \
        --tag $WRITETAG \
        --max-num-samples 1000 \
        --verbose
done 
