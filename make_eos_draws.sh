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
    python3 ../../draw_eos_uniform.py --num-draws $eos_per_dir 
    cd .. 
done
cd ..
