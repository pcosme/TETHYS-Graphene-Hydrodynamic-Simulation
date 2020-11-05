#!/bin/bash

DIRNAME=$(date +"DataFiles_BOUNDARY_LAYER_TRIAL_BATCH_FTCS_Slip_Length_1.5_with_Goduov_splitting_%F_%H:%M")
mkdir -p "$DIRNAME"

v=10
l=0.0
vis=0.0
cyc=0.0
for s in 20 40 60 
do 
	for vis in 0.05 0.1 0.2 0.6 1.0 1.2 1.5 
	do 
		echo "calculating S="$s
		echo "calculating Vf="$v
		echo "calculating visco="$vis
		./TETHYS_2D $s $v $l $vis $cyc 1
		mv -- *.dat "./$DIRNAME"
		mv -- *.h5 "./$DIRNAME"
	done
done 

mv -- *.log "./$DIRNAME"


