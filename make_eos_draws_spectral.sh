eos_per_dir=$1
num_mod_dirs=$2
eos_dir_name=$3

mkdir $eos_dir_name
cd $eos_dir_name
for index in $(seq -f '%06g' 1 $num_mod_dirs)
do
    dir=DRAWmod$eos_per_dir-$index
    mkdir $dir
    cd $dir
    # This is a bit of a problem, I don't really want to 
    # write code for python 2, but much of the code
    # I'm using is incompatible with python 3
    python3 ../../draw_eos_spectral.py --num-draws $eos_per_dir 
    cd .. 
done
cd ..
