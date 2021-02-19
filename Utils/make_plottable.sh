NUM_WEIGHTS=6
TAGS="parametric" "nonparametric"
for TAG in $TAGS
do
    python3 ~/parametric-eos-priors/Utils/correct_column_names.py \
        --posterior-file-path {$TAG}_post.csv \
        --output corrected_{$TAG}_post.csv \
        --num-weight-columns $NUM_WEIGHTS
done

condor_submit add_submit




