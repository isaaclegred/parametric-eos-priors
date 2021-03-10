#! /bin/bash
dir_index=$1
dirs_to_make=$2
eos_per_dir=$3

((min_index=dir_index*dirs_to_make))
eos_dir_name="gaussian_half_sigma_eos_draw_spectral"
((max_index=min_index + dirs_to_make - 1))

mkdir $eos_dir_name
cd $eos_dir_name
for raw_index in $(seq $min_index $max_index)
do
    # Should just be able to do seq -f %06g $min_index $max_index
    # but for some reason this doesn't work here
    index=$(printf "%06d" $raw_index)
    dir=DRAWmod$eos_per_dir-$index
    mkdir $dir
    cd $dir
    # This is a bit of a problem, I don't really want to 
    # write code for python 2, but much of the code
    # I'm using is incompatible with python 3
    python3 ../../draw_eos_spectral.py --num-draws $eos_per_dir --dir-index $index 
    cd .. 
done
cd ..
