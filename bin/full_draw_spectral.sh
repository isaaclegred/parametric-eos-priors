RUNDIR="$HOME/parametric-eos-priors/rundrawandsolve_spectralgauss/"
rm $RUNDIR/*
. make_draw_and_solve.sh "gaussian_eos_draw_spectral" "spectral gaussian" 100 10 200 $RUNDIR
