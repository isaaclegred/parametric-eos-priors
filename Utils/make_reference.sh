#!/bin/bash
MAIN_DIR=$(pwd)

PARAMS="Lambda(M=1.4)"
NEW_COVAR="R(M=1.4)"

TAGS="/corrected_nonparametric_post /corrected_parametric_post"

ADD_PARAMETERS=""
for PARAM in $PARAMS
do
    ADD_PARAMETERS=$ADD_PARAMETERS" --add-parameter "$PARAM 
done
echo $PARAMETER_LINES
for TAG in $TAGS
do
    POST_PATH=$MAIN_DIR/$TAG".csv"
    OUTPUT=$MAIN_DIR/$TAG"ROFLAMBDACCHYY.csv"
    python3 Utils/make_reference_csv.py --posterior-file-path $POST_PATH --output $OUTPUT --new-cov-name $NEW_COVAR $ADD_PARAMETERS
         
done
