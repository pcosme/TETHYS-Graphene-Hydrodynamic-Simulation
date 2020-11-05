#!/bin/bash

DIRNAME=$(date +"DataFiles_SINGLE_TEST_%F_%H:%M")
mkdir -p "$DIRNAME"
mkdir -p "$DIRNAME/extrema"
mkdir -p "$DIRNAME/electronics"

s=40
v=10
l=0.0
vis=0.0
cyc=0.0
echo "calculating S="$s
echo "calculating Vf="$v
echo "calculating Vc="$l
./TETHYS_1D $s $v $l $vis $cyc 1 

FILENAME1=$(find preview*.dat)
WORDCOUNT1=$(wc -l preview*.dat)
LINENUMBER1=${WORDCOUNT1% *}

./TimeSeries "$LINENUMBER1" "$FILENAME1" "$s" 

FILENAME2=$(find electro_*.dat)
WORDCOUNT2=$(wc -l electro_*.dat)
LINENUMBER2=${WORDCOUNT2% *}


./ElectronicAnalysis "$LINENUMBER2" "$FILENAME2" "$s" 

mv -- electro_S=*vF=*vis=*l=*.dat "./$DIRNAME/electronics"
mv -- extrema* "./$DIRNAME/extrema"
mv -- *.dat "./$DIRNAME"
mv -- *.h5 "./$DIRNAME"
mv -- *.log "./$DIRNAME"


