NUM_WEIGHTS=6
TAGS="parametric nonparametric"
for TAG in $TAGS
do
    # This is kinda silly, I think something needs to be refactored to make
    # everything make sense here, for example scripts should be different
    # from functions.
    python3 ~/parametric-eos-priors/Utils/correct_column_names.py \
        --posterior-file-path ${TAG}_post.csv \
        --output corrected_${TAG}_post.csv \
        --num-weight-columns $NUM_WEIGHTS
done
# I don't know about this either
cd Utils/Plotting
condor_submit add_submit
cd ../..



