//
// Created by pcosme on 04/02/2022.
//

#include "Cell1DLib.h"

CellHandler1D::CellHandler1D(int i, Fluid1D * ptr_fluid, StateVec * ptr_state) {
	index=i;
	fluid_ptr=ptr_fluid;
	U_ptr=ptr_state;
}

float CellHandler1D::VanLeer(int i) {
	float denom,numer,r=0,f=0;
	float limit=2.0f;
	float tolerance=1E-6;
	StateVec Numer = U_ptr[i]-U_ptr[i-1];
	StateVec Denom = U_ptr[i+1]-U_ptr[i];
	numer=Numer.v();
	denom=Denom.v();

	if(denom!=0){
		r=numer/denom;
		f=(r+abs(r))/(1+abs(r));
	} else{
		if(numer<=tolerance){
			f=0;
		} else{
			f=limit;
		}
	}
	return f;
}


float CellHandler1D::Roe(int i) {
	float denom,numer,r=0,f=0;
	float limit=1.0f;
	float tolerance=1E-5;
	StateVec Numer = U_ptr[i]-U_ptr[i-1];
	StateVec Denom = U_ptr[i+1]-U_ptr[i];
	numer=Numer.v();
	denom=Denom.v();

	if(denom!=0){
		r=numer/denom;
		f=max(0.0f,min(1.0f,r));
	} else{
		if(numer<=tolerance){
			f=0;
		} else{
			f=limit;
		}
	}
	return f;
}

StateVec CellHandler1D::TVD(char side, char edge) {
	StateVec Utvd(U_ptr[0]);
	switch(side) {
		case 'W' :
			switch(edge) {
				case 'L' :
					Utvd=U_ptr[index]   + 0.5f * VanLeer(index) * (U_ptr[index + 1] - U_ptr[index]);
					break;
				case 'R' :
					Utvd=U_ptr[index+1]   - 0.5f * VanLeer(index+1) * (U_ptr[index + 2] - U_ptr[index+1]);
					break;
				default: 0;
			}
			break;
		case 'E' :
			switch(edge) {
				case 'L' :
					Utvd=U_ptr[index-1]   + 0.5f * VanLeer(index-1) * (U_ptr[index] - U_ptr[index-1]);
					break;
				case 'R' :
					Utvd=U_ptr[index]   - 0.5f * VanLeer(index) * (U_ptr[index + 1] - U_ptr[index]);
					break;
				default: 0;
			}
			break;
		default: 0;
	}
	return Utvd;
}

StateVec CellHandler1D::UNO(char side, char edge) {
	StateVec Uuno(U_ptr[0]);
	return Uuno;
}

StateVec NumericalFlux::Average(Fluid1D *fluido, StateVec L, StateVec R) {
	StateVec Ureturn(fluido->DensityFlux(0.5f*(L+R)),fluido->VelocityFlux(0.5f*(L+R)));
	return Ureturn;
}

StateVec NumericalFlux::Central(Fluid1D *fluido, StateVec L, StateVec R) {
	StateVec Ureturn(fluido->DensityFlux(R)+fluido->DensityFlux(L),fluido->VelocityFlux(R)+fluido->VelocityFlux(L));
	float A =max(fluido->JacobianSpectralRadius(R),fluido->JacobianSpectralRadius(L));
	Ureturn.n() = Ureturn.n() - A*(R-L).n();
	Ureturn.v() = Ureturn.v() - A*(R-L).v();
	return 0.5f*Ureturn;
}

StateVec NumericalFlux::Characteristic(Fluid1D *fluido, StateVec L, StateVec R) {
	StateVec Ureturn(fluido->DensityFlux(R)+fluido->DensityFlux(L),fluido->VelocityFlux(R)+fluido->VelocityFlux(L));
	float u = 0.5f*(L+R).v();
	float A=(u/abs(u));
	Ureturn.n() = Ureturn.n() -A*(fluido->DensityFlux(R)-fluido->DensityFlux(L)  );
	Ureturn.v() = Ureturn.v() -A*(fluido->VelocityFlux(R)-fluido->VelocityFlux(L)  );
	return 0.5f*Ureturn;
}
