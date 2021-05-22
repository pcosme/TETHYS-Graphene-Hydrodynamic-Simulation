//
// Created by pcosme on 28/01/2021.
//


#include "includes/Grid2DLib.h"

GridPoint::GridPoint(int pos,int xpoints, int ypoints,bool mid) {
	IsMidGrid=mid;
	C=pos;
	size_x=xpoints;
	size_y=ypoints;
	if(IsMidGrid){
		stride=size_x-1;
	}else{
		stride=size_x;
	}
	FirstNeighbours();
	SecondNeighbours();
	ThirdNeighbours();
}

void GridPoint::FirstNeighbours() {
	if(IsMidGrid){
		stride=size_x-1;
		div_t divresult;
		divresult = div (C,stride);
		int j=divresult.quot;
		int i=divresult.rem;
		NE=i+1+(j+1)*(stride+1);
		NW=i+  (j+1)*(stride+1);
		SE=i+1    +j*(stride+1);
		SW=i      +j*(stride+1);
	}else{
		stride=size_x;
		div_t divresult;
		divresult = div (C,stride);
		int j=divresult.quot;
		int i=divresult.rem;
		NE=i      +j*(stride-1);
		NW=i-1    +j*(stride-1);
		SE=i  +(j-1)*(stride-1);
		SW=i-1+(j-1)*(stride-1);
	}
}


void GridPoint::SecondNeighbours() {
	if(IsMidGrid){
		stride=size_x-1;
	}else{
		stride=size_x;
	}
	N=C+stride;
	S=C-stride;
	E=C+1;
	W=C-1;
}

void GridPoint::ThirdNeighbours() {
	if(IsMidGrid){
		stride=size_x-1;
	}else{
		stride=size_x;
	}
	NE2=C+stride+1;
	NW2=C+stride-1;
	SE2=C-stride+1;
	SW2=C-stride-1;
}
