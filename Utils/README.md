To use these functions the first thing to do is to have posterior files as produced
from the eos-inf function `combine-obs`, I name them `parametric_post.csv` and
`nonparametric_post.csv` because I'm usually comparing nonparametric to parametric.

Put these in an empty directory somewhere, then copy this Utils folder into the
directory, as the same level as the posterior files.  Then edit the function
`Utils/Plotting/add-covariates` to point them to the correct draw directories.

Then call the function `Utils/make_plottable.sh` from the directory containing the
posterior files, this does two things.  It copies the posterior files and edits the
header to standardize it using the function `Utils/correct_column_names.py` (you'll
need to edit this function if your posterior files are named differently from
`parametric_post.csv, nonparametric_post.csv`, and then also change what functions
`Utils/Plotting/add-covariates` is looking for.  

Finally call the script `plot_corner.sh` which makes the corner plots, you shouldn't
have to edit anything, but in case something breaks, it is happening in the
script `Utils/Plotting/corner-covariates`, as there's also plotting settings there.

To review
(1) Follow eos-inf instruction all the way to combine-obs, do this for
a parametric and a nonparametric set of prior draws (you can do it with anything
but make sure to keep names consistent)

(2) Edit add-covariates to point to the correct prior draw directories

(3) Run make_plottable.sh

(4) Run plot_corner.sh


You can also make the pressure-density curves  by editing the file
`Utils/Plotting/test-quantiles`, in particular making sure again that you're
pointing to the correct draw directories, and either the name convention
with the post files is either kept or changed consistently. 
Then (again from the same directory as the
where the posterior files are) calling `Utils/plot_process.sh`, which will
submit the job to condor.  

