#!/bin/bash

dbdir="/home/isaac.legred/EOSDrawTest/eos_draw"
initeosnum=0
numeos=1100
eosperdir=50 

for eos in $(seq $initeosnum $(($initeosnum+$numeos-1)))
do
	echo $eos
	dirnum=$(($eos/$eosperdir))
	printf -v dirlabel "%06d" $dirnum
	printf -v eoslabel "%06d" $eos

	getnsprops eos-draw-$eoslabel.csv -v -p R,M,Lambda,I,Mb -m 3e6 -d $dbdir/DRAWmod$eosperdir-$dirlabel/ -o $dbdir/DRAWmod$eosperdir-$dirlabel/
	
	splitbranches macro-draw-$eoslabel.csv -v -d $dbdir/DRAWmod$eosperdir-$dirlabel/ -o $dbdir/DRAWmod$eosperdir-$dirlabel/ -f MACROdraw-$eoslabel -t "" 
	
	#plotprops macro-draw-$eoslabel.csv -v -p rhoc,R,Lambda,I,Mb -d $dbdir/DRAWmod$eosperdir-$dirlabel/ -o $dbdir/DRAWmod$eosperdir-$dirlabel/MACROdraw-$eoslabel/ -t draw-$eoslabel
done




