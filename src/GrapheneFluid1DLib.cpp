//
// Created by pcosme on 27/12/2020.
//

#include "GrapheneFluid1DLib.h"

GrapheneFluid1D::GrapheneFluid1D(SetUpParameters &input_parameters) : Fluid1D(input_parameters) {
	vel_fer = input_parameters.FermiVelocity;
	col_freq = input_parameters.CollisionFrequency;
	char buffer [50];
	sprintf (buffer, "S=%.2fvF=%.2fvis=%.2fl=%.2f", vel_snd, vel_fer, kin_vis, col_freq);
	file_infix = buffer;
}

GrapheneFluid1D::~GrapheneFluid1D(){
	delete Den;
	delete Vel ;
	delete Cur ;
	delete den_mid ;
	delete vel_mid ;
	delete DenCor ;
	delete VelCor ;
	delete CurCor ;
	delete vel_snd_arr ;
	delete GradVel ;
	delete grad_vel_mid ;
}

/*....................................................................*/

float GrapheneFluid1D::DensityFlux(float n,float v,float __attribute__((unused)) s){
	float f_1;
	f_1 = n * v;
	return f_1;
}

float GrapheneFluid1D::VelocityFlux(float n, float v, __attribute__((unused))float dv, float s){
	float f_2;
		f_2 = 0.25f * v * v + vel_fer * vel_fer * 0.5f * log(n) + 2.0f * s * s * sqrt(n);//- kin_vis * dv;
	return f_2;
}

float GrapheneFluid1D::DensitySource(__attribute__((unused)) float n, __attribute__((unused)) float v, __attribute__((unused)) float s){
	return 0.0f;
}

float GrapheneFluid1D::VelocitySource(__attribute__((unused)) float n, float v, __attribute__((unused)) float s){
	return -1.0f * col_freq * v ;
}

void GrapheneFluid1D::CflCondition(){
	dx = lengX / ( float ) ( Nx - 1 );
	float lambda;
	if(vel_snd<0.36f*vel_fer){
		lambda=1.2f*vel_fer;
	}else{
		lambda=1.97f*vel_snd + 0.5f*vel_fer;
	}
	dt = dx/lambda;
}