# To be run inside the eos_dir with incomplete directories
eos_per_dir=$1
num_mod_dirs=$2
eos_dir=$3
cd $eos_dir

for index in $(seq -f '%06g' 1 $num_mod_dirs)
do
    dir=DRAWmod$eos_per_dir-$index
    cd $dir
    if [ $(ls | wc -l) != $eos_per_dir ]
    then
        rm eos-draw-*.csv
        python3 ../../draw_eos_uniform.py --num-draws $eos_per_dir 
    fi 
    cd .. 
done
cd ..
