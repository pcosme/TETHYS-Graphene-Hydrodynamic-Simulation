#!/bin/bash

DIRNAME=$(date +"DataFiles_SINGLE_TEST_2D_%F_%H:%M")
mkdir -p "$DIRNAME"
mkdir -p "$DIRNAME/extrema"
mkdir -p "$DIRNAME/electronics"

s=40
v=10
l=0.0
vis=0.05
cyc=0.0
echo "calculating S="$s
echo "calculating Vf="$v
echo "calculating Vc="$l
./TETHYS_2D $s $v $l $vis $cyc 1

mv -- *.dat "./$DIRNAME"
mv -- *.h5 "./$DIRNAME"
mv -- *.log "./$DIRNAME"
