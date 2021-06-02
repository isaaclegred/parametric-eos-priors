#!/bin/bash

TAGS="corrected_all_miller"                                                                                             
MAIN_DIR=$(pwd)                                                                                      


for PRETAG in $TAGS
do TAG=$MAIN_DIR"/"$PRETAG"_post.csv"
   filter-samples \
       $TAG \
       $MAIN_DIR"/"$PRETAG"_post_allowed.csv" $MAIN_DIR"/"$PRETAG"_post_excluded".csv "data[:,cols.index('logweight_total')]>=-20"
done
