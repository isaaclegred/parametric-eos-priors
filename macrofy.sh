#!/bin/bash
# Trying to parallelize this
start_index=$1
# If we have access to a lot of procs we can make this smaller
# Total_dirs = num_procs * dirs_per_proc
dirs_to_make=$2


dbdir="/home/isaac.legred/parametric-eos-priors/production_eos_draw_spectral/"
eosperdir=200
(( total_eos_to_make=dirs_to_make*eosperdir ))
(( initeosnum=start_index*eosperdir*dirs_to_make ))



for eos in $(seq $initeosnum $(($initeosnum+$total_eos_to_make-1)))
do
	echo $eos
	dirnum=$(($eos/$eosperdir))
	printf -v dirlabel "%06d" $dirnum
	printf -v eoslabel "%06d" $eos

	getnsprops eos-draw-$eoslabel.csv -v -p R,M,Lambda,I,Mb -m 3e6 -d $dbdir/DRAWmod$eosperdir-$dirlabel/ -o $dbdir/DRAWmod$eosperdir-$dirlabel/
	
	splitbranches macro-draw-$eoslabel.csv -v -d $dbdir/DRAWmod$eosperdir-$dirlabel/ -o $dbdir/DRAWmod$eosperdir-$dirlabel/ -f MACROdraw-$eoslabel -t "" 
	
	#plotprops macro-draw-$eoslabel.csv -v -p rhoc,R,Lambda,I,Mb -d $dbdir/DRAWmod$eosperdir-$dirlabel/ -o $dbdir/DRAWmod$eosperdir-$dirlabel/MACROdraw-$eoslabel/ -t draw-$eoslabel
done




