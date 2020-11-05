#!/bin/bash


DIRNAME=$(date +"DataFiles_MAGNETIC_FIELD_TransverseCurrent_%F_%H:%M")
mkdir -p "$DIRNAME"

v=10
l=0.0
vis=0.0
cyc=0.0
for s in 10 20 40 60 
do 
	for cyc in 0.1 0.2 0.5 1.0 2.0 5.0 10.0
	do 
	echo "calculating S="$s
	echo "calculating Vf="$v
	echo "calculating cyclotron freq="$cyc
	./TETHYS_2D $s $v $l $vis $cyc 1
	mv -- *.dat "./$DIRNAME"
	mv -- *.h5 "./$DIRNAME"
	done
done 

mv -- *.log "./$DIRNAME"


