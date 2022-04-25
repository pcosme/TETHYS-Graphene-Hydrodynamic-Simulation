//
// Created by pcosme on 04/02/2022.
//

#include "Cell1DLib.h"

CellHandler1D::CellHandler1D(int i, int total, Fluid1D * ptr_fluid, StateVec * ptr_state) {
	index=i;
	size=total;
	fluid_ptr=ptr_fluid;
	U_ptr=ptr_state;
}

CellHandler1D::CellHandler1D(int i, int total, StateVec * ptr_state) {
	index=i;
	size=total;
	fluid_ptr= nullptr;
	U_ptr=ptr_state;
}

StateVec CellHandler1D::VanLeer(StateVec*Uin, int i) {
	StateVec Ureturn(Uin[0]);
	float denom,numer,r,f_vel=0,f_den=0;
	float tolerance=1E-6;

	StateVec Numer(Uin[i]);
	StateVec Denom(Uin[i]);
	if(i==0){
		Numer = Uin[0];
	} else{
		Numer = Uin[i]-Uin[i-1];
	}
	if(i==size-1){
		Denom = -1.0f*Uin[i];
	}else{
		Denom = Uin[i+1]-Uin[i];
	}


	/////////////////////////////VELOCITY
	numer=Numer.v();
	denom=Denom.v();
	r=numer/(denom+tolerance);
	f_vel=(r+abs(r))/(1+abs(r));

	/////////////////////////////Density
	numer=Numer.n();
	denom=Denom.n();
	r=numer/(denom+tolerance);

	f_den=(r+abs(r))/(1+abs(r));


	Ureturn.v()=f_vel;
	Ureturn.n()=f_den;
	return Ureturn;
}


StateVec CellHandler1D::Roe(StateVec*Uin, int i) {
	StateVec Ureturn(Uin[0]);
	float denom,numer,r,f_vel=0,f_den=0;
	float tolerance=1E-6;
	StateVec Numer = Uin[i]-Uin[i-1];
	StateVec Denom = Uin[i+1]-Uin[i];

	/////////////////////////////VELOCITY
	numer=Numer.v();
	denom=Denom.v();
	r=numer/(denom+tolerance);
	f_vel=max(0.0f,min(1.0f,r));
	/////////////////////////////Density
	numer=Numer.n();
	denom=Denom.n();
	r=numer/(denom+tolerance);
	f_den=max(0.0f,min(1.0f,r));

	Ureturn.v()=f_vel;
	Ureturn.n()=f_den;
	return Ureturn;
}


StateVec CellHandler1D::VanLeer(int i) {
	StateVec Ureturn(U_ptr[0]);
	float denom,numer,r,f_vel=0,f_den=0;
	float tolerance=1E-6;
	StateVec Numer = U_ptr[i]-U_ptr[i-1];
	StateVec Denom = U_ptr[i+1]-U_ptr[i];

	/////////////////////////////VELOCITY
	numer=Numer.v();
	denom=Denom.v();
	r=numer/(denom+tolerance);
	f_vel=(r+abs(r))/(1+abs(r));

	/////////////////////////////Density
	numer=Numer.n();
	denom=Denom.n();
	r=numer/(denom+tolerance);

	f_den=(r+abs(r))/(1+abs(r));


	Ureturn.v()=f_vel;
	Ureturn.n()=f_den;
	return Ureturn;
}


StateVec CellHandler1D::Roe(int i) {
	StateVec Ureturn(U_ptr[0]);
	float denom,numer,r,f_vel=0,f_den=0;
	float tolerance=1E-6;
	StateVec Numer = U_ptr[i]-U_ptr[i-1];
	StateVec Denom = U_ptr[i+1]-U_ptr[i];

	/////////////////////////////VELOCITY
	numer=Numer.v();
	denom=Denom.v();
	r=numer/(denom+tolerance);
	f_vel=max(0.0f,min(1.0f,r));
	/////////////////////////////Density
	numer=Numer.n();
	denom=Denom.n();
	r=numer/(denom+tolerance);
	f_den=max(0.0f,min(1.0f,r));

	Ureturn.v()=f_vel;
	Ureturn.n()=f_den;
	return Ureturn;
}

StateVec CellHandler1D::TVD(StateVec * Uin,int pos,char side, char edge) {
	StateVec Utvd(Uin[pos]);
	switch(side) {
		case 'E' :
			switch(edge) {
				case 'L' :
					Utvd=Uin[pos]   + 0.5f * VanLeer(Uin, pos) * (Uin[pos + 1] - Uin[pos]); //TODO passar a usar o limitador que só depende da posicao
					break;
				case 'R' :
					Utvd=Uin[pos+1]   - 0.5f * VanLeer(Uin, pos + 1) * (Uin[pos + 2] - Uin[pos + 1]);
					break;
				default: 0;
			}
			break;
		case 'W' :
			switch(edge) {
				case 'L' :
					Utvd=Uin[pos-1]   + 0.5f * VanLeer(Uin, pos - 1) * (Uin[pos] - Uin[pos - 1]);
					break;
				case 'R' :
					Utvd=Uin[pos]   - 0.5f * VanLeer(Uin, pos) * (Uin[pos + 1] - Uin[pos]);
					break;
				default: 0;
			}
			break;
		default: 0;
	}
	return Utvd;
}


StateVec CellHandler1D::TVD(char side, char edge) {
	StateVec Utvd(U_ptr[index]);
	int pos=index;
	switch(side) {
		case 'E' :
			switch(edge) {
				case 'L' :
					Utvd=U_ptr[pos]   + 0.5f * VanLeer(U_ptr, pos) * (U_ptr[pos + 1] - U_ptr[pos]);
					break;
				case 'R' :   //problematico em index==Nx-2
					if(pos==size-2){
						Utvd = U_ptr[pos + 1] - 0.5f * VanLeer(U_ptr, pos + 1) * ( -1.0f* U_ptr[pos + 1]);
					}else {
						Utvd = U_ptr[pos + 1] - 0.5f * VanLeer(U_ptr, pos + 1) * (U_ptr[pos + 2] - U_ptr[pos + 1]);
					}
					break;
				default: 0;
			}
			break;
		case 'W' :
			switch(edge) {
				case 'L' :
					Utvd=U_ptr[pos-1]   + 0.5f * VanLeer(U_ptr, pos - 1) * (U_ptr[pos] - U_ptr[pos - 1]);
					break;
				case 'R' :
					Utvd=U_ptr[pos]   - 0.5f * VanLeer(U_ptr, pos) * (U_ptr[pos + 1] - U_ptr[pos]);
					break;
				default: 0;
			}
			break;
		default: 0;
	}
	return Utvd;
}


StateVec NumericalFlux::Average(Fluid1D *fluido, StateVec L, StateVec R) {
	StateVec Ureturn{};
	Ureturn.n()=fluido->DensityFlux(0.5f*(L+R));
	Ureturn.v()=fluido->VelocityFlux(0.5f*(L+R));
	//Ureturn.S()=fluido->GetVelSnd();
	return Ureturn;
}

StateVec NumericalFlux::Central(Fluid1D *fluido, StateVec L, StateVec R) {
	StateVec Ureturn(L);
	Ureturn.n()=fluido->DensityFlux(R)+fluido->DensityFlux(L);
	Ureturn.v()=fluido->VelocityFlux(R)+fluido->VelocityFlux(L);
	//Ureturn.S()=fluido->GetVelSnd();
	float A =max(fluido->JacobianSpectralRadius(R),fluido->JacobianSpectralRadius(L));
	Ureturn.n() = Ureturn.n() - A*(R-L).n();
	Ureturn.v() = Ureturn.v() - A*(R-L).v();
	return 0.5f*Ureturn;
}

StateVec NumericalFlux::Characteristic(Fluid1D *fluido,  StateVec L, StateVec R) {
	StateVec Ureturn(fluido->DensityFlux(R)+fluido->DensityFlux(L),fluido->VelocityFlux(R)+fluido->VelocityFlux(L));

	StateVec u{};
	u= 0.5f*(L+R);

	float Fden,Fvel;

	Fden = fluido->DensityFlux(R)-fluido->DensityFlux(L);
	Fvel = fluido->VelocityFlux(R)-fluido->VelocityFlux(L);

	Ureturn.n() = Ureturn.n() - ( fluido->JacobianSignum(u,"11")*Fden + fluido->JacobianSignum(u,"12")*Fvel );
	Ureturn.v() = Ureturn.v() - ( fluido->JacobianSignum(u,"21")*Fden + fluido->JacobianSignum(u,"22")*Fvel );

	return 0.5f*Ureturn;
}

