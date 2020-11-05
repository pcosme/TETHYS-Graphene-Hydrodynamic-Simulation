#include "BoundaryLib.h"

#ifndef MAT_PI
#    define MAT_PI 3.14159265358979323846
#endif

#ifndef MAT_EULER
#    define MAT_EULER 2.71828182845905
#endif


using namespace H5;
using namespace std;


void BoundaryCondition::XFree(Fluid1D& fluid_class){
	int nx=fluid_class.SizeX();
	fluid_class.Den[0] = fluid_class.Den[1];
	fluid_class.Den[nx - 1] =  fluid_class.Den[nx - 2];
	fluid_class.Vel[0] = fluid_class.Vel[1];
	fluid_class.Vel[nx - 1] = fluid_class.Vel[nx - 2];
}
void BoundaryCondition::XPeriodic(Fluid1D& fluid_class){
	int nx=fluid_class.SizeX();
	fluid_class.Den[0] = fluid_class.Den[nx - 2];
	fluid_class.Den[nx - 1] = fluid_class.Den[1];
	fluid_class.Vel[0] = fluid_class.Vel[nx - 2];
	fluid_class.Vel[nx - 1] = fluid_class.Vel[1];
}
void BoundaryCondition::XFree(Fluid2D& fluid_class){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	
	for(int j=0; j < ny; j++){
		int left= 0 + j * nx;
		int right= nx - 1 + j * nx;
		fluid_class.Den[left]=fluid_class.Den[left + 1];
		fluid_class.Den[right]=fluid_class.Den[right - 1];			//free density at x=L
		fluid_class.FlxY[left] = fluid_class.FlxY[left + 1]; 					//flux only on x at x=0
		fluid_class.FlxY[right] =  fluid_class.FlxY[right - 1] ;					//idem at x=L
		fluid_class.FlxX[left] = fluid_class.FlxX[left + 1] ;			//free flux at x=0
		fluid_class.FlxX[right] =  fluid_class.FlxX[right - 1];

		//TODO REVER ISTO ENTAO NAO TA BEM ABERTO
		/*
		fluid_class.FlxY[left] = 0.0f; 					//flux only on x at x=0
		fluid_class.FlxY[right] = 0.0f ;					//idem at x=L
		fluid_class.FlxX[left] = fluid_class.FlxX[left + 1] * pow(fluid_class.Den[left + 1], -1.5f);			//free flux at x=0
		fluid_class.FlxX[right] =  fluid_class.FlxX[right - 1];*/
	}	
}
void BoundaryCondition::XFreeLeft(Fluid2D &fluid_class) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for(int j=0; j < ny; j++){
		int left= 0 + j * nx;
		fluid_class.Den[left]=fluid_class.Den[left + 1];
		fluid_class.FlxY[left] = 0.0f; 					//flux only on x at x=0
		fluid_class.FlxX[left] = fluid_class.FlxX[left + 1] * pow(fluid_class.Den[left + 1], -1.5f);			//free flux at x=0
	}
}
void BoundaryCondition::XFreeRight(Fluid2D &fluid_class) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for(int j=0; j < ny; j++){
		int right= nx - 1 + j * nx;
		fluid_class.Den[right]=fluid_class.Den[right - 1];			//free density at x=L
		fluid_class.FlxY[right] = 0.0f ;					//idem at x=L
		fluid_class.FlxX[right] =  fluid_class.FlxX[right - 1];
	}
}


void BoundaryCondition::XPeriodic(Fluid2D& fluid_class){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	
	for(int j=0; j < ny; j++){
		int left= 0 + j * nx;
		int right= nx - 1 + j * nx;
		fluid_class.Den[left]=fluid_class.Den[right - 1];
		fluid_class.Den[right]=fluid_class.Den[1 + j * nx];
		fluid_class.FlxY[left] = 0.0; 					//flux only on x at x=0
		fluid_class.FlxY[right] = 0.0 ;					//idem at x=L
		fluid_class.FlxX[left] = fluid_class.FlxX[right - 1];
		fluid_class.FlxX[right] =  fluid_class.FlxX[left + 1];
	}	
}
void BoundaryCondition::YFree(Fluid2D& fluid_class){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int i=0; i < nx; i++){
		int bottom=i; //i+0*nx
		int top= i + (ny - 1) * nx;
		fluid_class.Den[bottom] = fluid_class.Den[bottom + nx];
		fluid_class.FlxX[bottom] = fluid_class.FlxX[bottom + nx];
		fluid_class.FlxY[bottom] = fluid_class.FlxY[bottom + nx];
		fluid_class.Den[top] = fluid_class.Den[top - nx];
		fluid_class.FlxX[top] = fluid_class.FlxX[top - nx];
		fluid_class.FlxY[top] = fluid_class.FlxY[top - nx];

	}	 	
}
void BoundaryCondition::YFreeTop(Fluid2D &fluid_class) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int i=0; i < nx; i++){
		int top= i + (ny - 1) * nx;
		fluid_class.Den[top] = fluid_class.Den[top - nx];
		fluid_class.FlxX[top] = fluid_class.FlxX[top - nx];
		fluid_class.FlxY[top] = fluid_class.FlxY[top - nx];
	}
}
void BoundaryCondition::YFreeBottom(Fluid2D &fluid_class) {
	int nx=fluid_class.SizeX();
	for (int i=0; i < nx; i++){
		int bottom=i; //i+0*nx
		fluid_class.Den[bottom] = fluid_class.Den[bottom + nx];
		fluid_class.FlxX[bottom] = fluid_class.FlxX[bottom + nx];
		fluid_class.FlxY[bottom] = fluid_class.FlxY[bottom + nx];
	}
}
void BoundaryCondition::YPeriodic(Fluid2D& fluid_class){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int i=0; i < nx; i++){
		int bottom=i; //i+0*nx
		int top= i + (ny - 1) * nx;
		fluid_class.Den[bottom] = fluid_class.Den[top - nx];
		fluid_class.FlxX[bottom] = fluid_class.FlxX[top - nx];
		fluid_class.FlxY[bottom] = fluid_class.FlxY[top - nx];
		fluid_class.Den[top] = fluid_class.Den[bottom + nx];
		fluid_class.FlxX[top] = fluid_class.FlxX[bottom + nx];
		fluid_class.FlxY[top] = fluid_class.FlxY[bottom + nx];
	}
}

void BoundaryCondition::YClosedFreeSlip(Fluid2D& fluid_class){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int i=0; i < nx; i++){
		int bottom=i; //i+0*nx
		int top= i + (ny - 1) * nx;
		fluid_class.Den[bottom] = fluid_class.Den[bottom + nx];
		fluid_class.FlxX[bottom] = fluid_class.FlxX[top - nx];
		fluid_class.FlxY[bottom] = 0.0f;
		fluid_class.Den[top] = fluid_class.Den[top - nx];
		fluid_class.FlxX[top] = fluid_class.FlxX[bottom + nx];
		fluid_class.FlxY[top] =  0.0f;
	}
}
void BoundaryCondition::YClosedNoSlip(Fluid2D& fluid_class){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int i=0; i < nx; i++){
		int bottom=i; //i+0*nx
		int top= i + (ny - 1) * nx;
		fluid_class.Den[bottom] = fluid_class.Den[bottom + nx];
		fluid_class.FlxX[bottom] =  0.0f;
		fluid_class.FlxY[bottom] =  0.0f;
		fluid_class.Den[top] = fluid_class.Den[top - nx];
		fluid_class.FlxX[top] =  0.0f;
		fluid_class.FlxY[top] =  0.0f;
	}
}




void DirichletBoundaryCondition::Density(Fluid1D& fluid_class, float left, float right){
	int nx=fluid_class.SizeX();
	fluid_class.Den[0] = left;
	fluid_class.Den[nx - 1] = right;
}
void DirichletBoundaryCondition::VelocityX(Fluid1D& fluid_class, float left, float right){
	int nx=fluid_class.SizeX();
	fluid_class.Vel[0] = left;
	fluid_class.Vel[nx - 1] = right;
}  
void DirichletBoundaryCondition::Density(Fluid2D& fluid_class, float left, float right, float top, float bottom){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		fluid_class.Den[j * nx] = left;
		fluid_class.Den[nx - 1 + j * nx] = right;
	}
	for (int i=0; i < nx; i++){
		fluid_class.Den[i + (ny - 1) * nx] = top;
		fluid_class.Den[i] = bottom;
	}
} 

void DirichletBoundaryCondition::MassFluxX(Fluid2D& fluid_class, float left, float right, float top, float bottom){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		fluid_class.FlxX[j * nx] = left;
		fluid_class.FlxX[nx - 1 + j * nx] = right;
	}
	for (int i=0; i < nx; i++){
		fluid_class.FlxX[i + (ny - 1) * nx] = top;
		fluid_class.FlxX[i] = bottom;
	}
} 
void DirichletBoundaryCondition::MassFluxY(Fluid2D& fluid_class, float left, float right, float top, float bottom){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		fluid_class.FlxY[j * nx] = left;
		fluid_class.FlxY[nx - 1 + j * nx] = right;
	}
	for (int i=0; i < nx; i++){
		fluid_class.FlxY[i + (ny - 1) * nx] = top;
		fluid_class.FlxY[i ] = bottom;
	}
}

void DirichletBoundaryCondition::Jet(Fluid2D &fluid_class, float left, float left_width, float right, float right_width) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	int n_width_left= static_cast<int>(static_cast<float>(ny) * left_width);
	int n_width_right= static_cast<int>(static_cast<float>(ny) * right_width);
	for (int j=0; j < ny; j++){
		if( j>=(ny-n_width_left)/2 && j<= (ny+n_width_left)/2){
			fluid_class.FlxX[0 + j * nx] = left;
		} else{
			fluid_class.FlxX[0 + j * nx] = 0.0f;
		}
		if( j>=(ny-n_width_right)/2 && j<= (ny+n_width_right)/2){
			fluid_class.FlxX[nx - 1 + j * nx] = right;
		} else{
			fluid_class.FlxX[nx - 1 + j * nx] = 0.0f;
		}
	}
}

void DirichletBoundaryCondition::DensityRight(Fluid2D &fluid_class, float right) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		fluid_class.Den[nx - 1 + j * nx] = right;
	}
}

void DirichletBoundaryCondition::MassFluxXRight(Fluid2D &fluid_class, float right) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		fluid_class.FlxX[nx - 1 + j * nx] = right;
	}
}

void DirichletBoundaryCondition::MassFluxYRight(Fluid2D &fluid_class, float right) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		fluid_class.FlxY[nx - 1 + j * nx] = right;
	}
}

void DirichletBoundaryCondition::DensityLeft(Fluid2D &fluid_class, float left) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		fluid_class.Den[j * nx] = left;
	}
}

void DirichletBoundaryCondition::MassFluxXLeft(Fluid2D &fluid_class, float left) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		fluid_class.FlxX[j * nx] = left;
	}
}

void DirichletBoundaryCondition::MassFluxYLeft(Fluid2D &fluid_class, float left) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		fluid_class.FlxY[j * nx] = left;
	}
}

void DirichletBoundaryCondition::DensityTop(Fluid2D &fluid_class, float top) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int i=0; i < nx; i++){
		fluid_class.Den[i + (ny - 1) * nx] = top;
	}
}

void DirichletBoundaryCondition::MassFluxXTop(Fluid2D &fluid_class, float top) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int i=0; i < nx; i++){
		fluid_class.FlxX[i + (ny - 1) * nx] = top;
	}
}

void DirichletBoundaryCondition::MassFluxYTop(Fluid2D &fluid_class, float top) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int i=0; i < nx; i++){
		fluid_class.FlxY[i + (ny - 1) * nx] = top;
	}
}

void DirichletBoundaryCondition::DensityBottom(Fluid2D &fluid_class, float bottom) {
	int nx=fluid_class.SizeX();
	for (int i=0; i < nx; i++){
		fluid_class.Den[i] = bottom;
	}
}

void DirichletBoundaryCondition::MassFluxXBottom(Fluid2D &fluid_class, float bottom) {
	int nx=fluid_class.SizeX();
	for (int i=0; i < nx; i++){
		fluid_class.FlxX[i] = bottom;
	}
}



void DirichletBoundaryCondition::MassFluxYBottom(Fluid2D &fluid_class, float bottom) {
	int nx=fluid_class.SizeX();
	for (int i=0; i < nx; i++){
		fluid_class.FlxY[i] = bottom;
	}
}

void DyakonovShurBoundaryCondition::DyakonovShurBc(GrapheneFluid1D& fluid_class) {
	int nx=fluid_class.SizeX();
	fluid_class.Den[0] = 1.0f;
	fluid_class.Vel[0] = fluid_class.Vel[1];
	fluid_class.Den[nx - 1] = fluid_class.Den[nx - 2];
	fluid_class.Vel[nx - 1] = 1.0f / fluid_class.Den[nx - 1];
}

void DyakonovShurBoundaryCondition::DyakonovShurBc(GrapheneFluid2D& fluid_class) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	 
	for(int j=0; j < ny; j++){
		fluid_class.Den[0 + j * nx]=1.0f;			//constant density at x=0
		fluid_class.Den[nx - 1 + j * nx]=fluid_class.Den[nx - 2 + j * nx]; 			//free density at x=L
		fluid_class.FlxX[0 + j * nx] = fluid_class.FlxX[1 + j * nx] * pow(fluid_class.Den[1 + j * nx], -1.5f);			//free flux at x=0
		fluid_class.FlxX[nx - 1 + j * nx] = sqrt(fluid_class.Den[nx - 1 + j * nx]);	//constant current at x=L (flux equals mass)
		fluid_class.FlxY[0 + j * nx] = 0.0f; 					//flux only on x at x=0
		fluid_class.FlxY[nx - 1 + j * nx] = 0.0f ;					//idem at x=L
	}	
}


void RobinBoundaryCondition::SlipLength(Fluid2D &fluid_class, float slip_length) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	float dy=fluid_class.GetDy();
	for (int i=0; i < nx; i++){
		float aux_top,aux_bot, dn_top,dn_bot,l_top,l_bot;
		int bottom=i; //i+0*nx
		int top= i + (ny - 1) * nx;
		dn_bot = (-3.0f*fluid_class.Den[i+0*nx]+4.0f*fluid_class.Den[i+1*nx]-1.0f*fluid_class.Den[i+2*nx])/(2.0f*dy);
		dn_top = (-3.0f*fluid_class.Den[i+(ny - 1)*nx]+4.0f*fluid_class.Den[i+(ny - 2)*nx]-1.0f*fluid_class.Den[i+(ny - 3)*nx])/(2.0f*dy);
		dn_bot = dn_bot/sqrt(fluid_class.Den[bottom]);
		dn_top = dn_top/sqrt(fluid_class.Den[top]);
		l_top = slip_length/(1.0f+slip_length*dn_top);
		l_bot = slip_length/(1.0f-slip_length*dn_bot);
		aux_top = l_top/(2.0f*dy+3.0f*l_top);
		aux_bot = l_bot/(2.0f*dy+3.0f*l_bot);
		fluid_class.Den[bottom] = fluid_class.Den[bottom + nx];
		fluid_class.Den[top] = fluid_class.Den[top - nx];
		fluid_class.FlxY[top] =  0.0f;
		fluid_class.FlxY[bottom] = 0.0f;
		fluid_class.FlxX[bottom] = aux_bot*(4.0f*fluid_class.FlxX[i+1*nx]-1.0f*fluid_class.FlxX[i+2*nx]);
		fluid_class.FlxX[top] = aux_top*(4.0f*fluid_class.FlxX[i+(ny-2)*nx]-1.0f*fluid_class.FlxX[i+(ny-3)*nx]);
	}
}

void RobinBoundaryCondition::SlipLengthTop(Fluid2D &fluid_class, float slip_length) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	float dy=fluid_class.GetDy();
	for (int i=0; i < nx; i++){
		float aux_top, dn_top,l_top;
		int top= i + (ny - 1) * nx;
		dn_top = (-3.0f*fluid_class.Den[i+(ny - 1)*nx]+4.0f*fluid_class.Den[i+(ny - 2)*nx]-1.0f*fluid_class.Den[i+(ny - 3)*nx])/(2.0f*dy);
		dn_top = dn_top/sqrt(fluid_class.Den[top]);
		l_top = slip_length/(1.0f+slip_length*dn_top);
		aux_top = l_top/(2.0f*dy+3.0f*l_top);
		fluid_class.Den[top] = fluid_class.Den[top - nx];
		fluid_class.FlxY[top] =  0.0f;
		fluid_class.FlxX[top] = aux_top*(4.0f*fluid_class.FlxX[i+(ny-2)*nx]-1.0f*fluid_class.FlxX[i+(ny-3)*nx]);
	}
}

void RobinBoundaryCondition::SlipLengthBottom(Fluid2D &fluid_class, float slip_length) {
	int nx=fluid_class.SizeX();
	float dy=fluid_class.GetDy();
	for (int i=0; i < nx; i++){
		float aux_bot, dn_bot, l_bot;
		int bottom=i; //i+0*nx
		dn_bot = (-3.0f*fluid_class.Den[i+0*nx]+4.0f*fluid_class.Den[i+1*nx]-1.0f*fluid_class.Den[i+2*nx])/(2.0f*dy);
		dn_bot = dn_bot/sqrt(fluid_class.Den[bottom]);
		l_bot = slip_length/(1.0f-slip_length*dn_bot);
		aux_bot = l_bot/(2.0f*dy+3.0f*l_bot);
		fluid_class.Den[bottom] = fluid_class.Den[bottom + nx];
		fluid_class.FlxY[bottom] = 0.0f;
		fluid_class.FlxX[bottom] = aux_bot*(4.0f*fluid_class.FlxX[i+1*nx]-1.0f*fluid_class.FlxX[i+2*nx]);
	}
}
