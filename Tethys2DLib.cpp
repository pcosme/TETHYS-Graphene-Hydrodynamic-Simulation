// 2D version

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <ctime>

#include <string>




#include "TethysLib.h"
#include "Tethys2DLib.h"
#include <H5Cpp.h>

using namespace H5;
using namespace std;


#ifndef MAT_PI
#	define MAT_PI 3.14159265358979323846
#endif

Fluid2D::Fluid2D(int size_nx, int size_ny, float sound_velocity, float shear_viscosity) : TETHYSBase{size_nx, size_ny, 2}{
	Nx = size_nx;
	Ny = size_ny;
	vel_snd =sound_velocity;
	kin_vis =shear_viscosity;
	char buffer [50];
	sprintf (buffer, "south=%.2fvis=%.2f", vel_snd, kin_vis);
	file_infix = buffer;
	// main grid variables Nx*Ny
	Den 		= new float[Nx * Ny]();
	VelX 		= new float[Nx * Ny]();
	VelY 		= new float[Nx * Ny]();
	FlxX 		= new float[Nx * Ny]();
	FlxY 		= new float[Nx * Ny]();
	CurX 		= new float[Nx * Ny]();
	CurY 		= new float[Nx * Ny]();
	vel_snd_arr	= new float[Nx*Ny]();

	lap_flxX = new float[Nx*Ny](); //new grids for the laplacians
	lap_flxY = new float[Nx*Ny](); //in fact they could be smaller but thiw way they are just 0 at the borders who do not evolve

	// 1st Aux. Grid variables (Nx-1)*(Ny-1)
	den_mid		= new float[(Nx-1)*(Ny-1)]();  
	flxX_mid	= new float[(Nx-1)*(Ny-1)]();
	flxY_mid	= new float[(Nx-1)*(Ny-1)]();

}
	
Fluid2D::~Fluid2D(){
	delete [] Den;
	delete [] VelX;
	delete [] VelY;
	delete [] FlxX;
	delete [] FlxY;
	delete [] CurX;
	delete [] CurY;
	delete [] den_mid;
	delete [] flxX_mid;
	delete [] flxY_mid;
	delete [] lap_flxX;
	delete [] lap_flxY;
	delete [] vel_snd_arr;
}



void Fluid2D::SetSound(){ 
	for(int i = 0; i<Nx  ;i++){
		for(int j=0; j<Ny ; j++){
			//vel_snd_arr[i+j*Nx]=Sound_Velocity_Anisotropy(i,dx,j,dy, vel_snd);
			vel_snd_arr[i+j*Nx]=vel_snd;
		}
	}
}


float Fluid2D::GetVelSnd(){ return vel_snd; }
void Fluid2D::SetVelSnd(float x){ vel_snd=x; }
float Fluid2D::GetKinVis(){ return kin_vis; }
void Fluid2D::SetKinVis(float x){ kin_vis=x;}
float Fluid2D::GetDx(){return dx;}
void Fluid2D::SetDx(float x){ dx=x;} 
float Fluid2D::GetDy(){return dy;}
void Fluid2D::SetDy(float x){ dy=x;}
float Fluid2D::GetDt(){return dt;}
void Fluid2D::SetDt(float x){ dt=x;}

void Fluid2D::InitialCondRand(){
	srand (static_cast<unsigned int>(time(NULL)));
	for (int i = 0; i < Nx; i++ ){
		for (int j=0; j<Ny; j++){
		float noise = (float) rand()/ (float) RAND_MAX ;
			Den[i + j * Nx] = 1.0f + 0.005f * (noise - 0.5f);
		}
	}
}

void Fluid2D::InitialCondTest(){
	for (int i = 0; i < Nx; i++ ){
		for (int j=0; j<Ny; j++){
			float densi;
			if(i>=80&&i<=120&&j>=80&&j<=120){
			densi=0.2f;
			}
			else{
			densi=0.0f;
			}
			Den[i + j * Nx] = 1.0f + densi;
			VelX[i + j * Nx] = 0.1f;
		}
	}
}


void Fluid2D::MassFluxToVelocity(){
	for(int c=0; c <= Nx * Ny - 1; c++){
		VelX[c]= FlxX[c] / Den[c];
		VelY[c]= FlxY[c] / Den[c];
		CurX[c] = VelX[c] * Den[c];
		CurY[c] = VelY[c] * Den[c];
	}
}

void Fluid2D::Richtmyer(){
	// TODO throw exception to cath NAN or INF values
		int northeast,northwest,southeast,southwest;
		float den_north, den_south ,den_east ,den_west, px_north, px_south, px_east, px_west, py_north, py_south, py_east, py_west,m_east,m_west,m_north,m_south;
		//k=i+j*Nx
		for(int ks=0; ks<=Nx*Ny-Nx-Ny; ks++){ //correr todos os pontos da grelha secundaria de den_mid
			div_t divresult;
			divresult = div (ks,Nx-1);
			int j=divresult.quot;
			int i=divresult.rem;

			northeast=i+1+(j+1)*Nx; //mal  ->i+1,j+1 Prin  kPrin = i+j*Nx
			northwest=i+(j+1)*Nx;   //mal  ->i,j+1   Prin
			southeast=i+1+j*Nx;    //mal  ->i+1,j   Prin
			southwest=i+j*Nx;      //mal  ->i,j     Prin
		
				den_north = 0.5f*(Den[northeast] + Den[northwest]);
				den_south = 0.5f*(Den[southeast] + Den[southwest]);
				den_east = 0.5f*(Den[northeast] + Den[southeast]);
				den_west = 0.5f*(Den[northwest] + Den[southwest]);

				px_north = 0.5f*(FlxX[northeast] + FlxX[northwest]);
				px_south = 0.5f*(FlxX[southeast] + FlxX[southwest]);
				px_east = 0.5f*(FlxX[northeast] + FlxX[southeast]);
				px_west = 0.5f*(FlxX[northwest] + FlxX[southwest]);
				
				py_north = 0.5f*(FlxY[northeast] + FlxY[northwest]);
				py_south = 0.5f*(FlxY[southeast] + FlxY[southwest]);
				py_east = 0.5f*(FlxY[northeast] + FlxY[southeast]);
				py_west = 0.5f*(FlxY[northwest] + FlxY[southwest]);

				//posso definir aqui a "massa" nos 4 ponto s cardeais
				m_east=pow(den_east,1.5f); // e assim sucessivamente m_west m_north m_south que depois sao reutilizadeas nos 12 fluxos
				m_west=pow(den_west,1.5f);
				m_north=pow(den_north,1.5f);
				m_south=pow(den_south,1.5f);
				den_mid[ks] = 0.25f*(Den[southwest] + Den[southeast] + Den[northwest] + Den[northeast]) // How shall we include vel_snd_arr ?
								-0.5f*(dt/dx)*(
									DensityFluxX(den_east, px_east, py_east,m_east,vel_snd)-
									DensityFluxX(den_west, px_west, py_west,m_west,vel_snd))
									-0.5f*(dt/dy)*(
									DensityFluxY(den_north, px_north, py_north,m_north,vel_snd)-
									DensityFluxY(den_south, px_south, py_south,m_south,vel_snd));
				flxX_mid[ks] = 0.25f*(FlxX[southwest] + FlxX[southeast] + FlxX[northwest] + FlxX[northeast])
								-0.5f*(dt/dx)*(
									MassFluxXFluxX(den_east, px_east, py_east,m_east,vel_snd)-
									MassFluxXFluxX(den_west, px_west, py_west,m_west,vel_snd))
								-0.5f*(dt/dy)*(
									MassFluxXFluxY(den_north, px_north, py_north,m_north,vel_snd)-
									MassFluxXFluxY(den_south, px_south, py_south,m_south,vel_snd));
				flxY_mid[ks] = 0.25f*(FlxY[southwest] + FlxY[southeast] + FlxY[northwest] + FlxY[northeast])
								-0.5f*(dt/dx)*(
									MassFluxYFluxX(den_east, px_east, py_east,m_east,vel_snd)-
									MassFluxYFluxX(den_west, px_west, py_west,m_west,vel_snd))
								-0.5f*(dt/dy)*(
									MassFluxYFluxY(den_north, px_north, py_north,m_north,vel_snd)-
									MassFluxYFluxY(den_south, px_south, py_south,m_south,vel_snd));
	//		}						
		}
		for(int kp=1+Nx; kp<=Nx*Ny-Nx-2; kp++){ //correr a grelha principal evitando as fronteiras
			div_t divresult;
			divresult = div (kp,Nx);
			int j=divresult.quot;
			int i=divresult.rem;


			if( kp%Nx!=Nx-1 && kp%Nx!=0){

				northeast=i+j*(Nx-1);     //mal  ->i,j Sec  kSec = i+j*(Nx-1)
				northwest=i-1+j*(Nx-1);   //mal  ->i-1,j Sec
				southeast=i+(j-1)*(Nx-1);  //mal  ->i,j-1 Sec
				southwest=i-1+(j-1)*(Nx-1); //mal ->i-1,j-1 Sec
				
				den_north = 0.5f*(den_mid[northeast]+den_mid[northwest]);
				den_south = 0.5f*(den_mid[southeast]+den_mid[southwest]);
				den_east = 0.5f*(den_mid[northeast]+den_mid[southeast]);
				den_west = 0.5f*(den_mid[northwest]+den_mid[southwest]);

				px_north = 0.5f*(flxX_mid[northeast]+flxX_mid[northwest]);
				px_south = 0.5f*(flxX_mid[southeast]+flxX_mid[southwest]);
				px_east = 0.5f*(flxX_mid[northeast]+flxX_mid[southeast]);
				px_west = 0.5f*(flxX_mid[northwest]+flxX_mid[southwest]);
				
				py_north = 0.5f*(flxY_mid[northeast]+flxY_mid[northwest]);
				py_south = 0.5f*(flxY_mid[southeast]+flxY_mid[southwest]);
				py_east = 0.5f*(flxY_mid[northeast]+flxY_mid[southeast]);
				py_west = 0.5f*(flxY_mid[northwest]+flxY_mid[southwest]);

				//posso definir aqui a "massa" nos 4 ponto s cardeais
				m_east=pow(den_east,1.5f); // e assim sucessivamente m_west m_north m_south que depois sao reutilizadeas nos 12 fluxos
				m_west=pow(den_west,1.5f);
				m_north=pow(den_north,1.5f);
				m_south=pow(den_south,1.5f);

				Den[kp] = Den[kp]-(dt/dx)*(
									DensityFluxX(den_east, px_east, py_east,m_east,vel_snd)-
									DensityFluxX(den_west, px_west, py_west,m_west,vel_snd))
									-(dt/dy)*(
									DensityFluxY(den_north, px_north, py_north,m_north,vel_snd)-
									DensityFluxY(den_south, px_south, py_south,m_south,vel_snd));
				FlxX[kp] = FlxX[kp]-(dt/dx)*(
									MassFluxXFluxX(den_east, px_east, py_east,m_east,vel_snd)-
									MassFluxXFluxX(den_west, px_west, py_west,m_west,vel_snd))
									-(dt/dy)*(
									MassFluxXFluxY(den_north, px_north, py_north,m_north,vel_snd)-
									MassFluxXFluxY(den_south, px_south, py_south,m_south,vel_snd));
				FlxY[kp] = FlxY[kp]-(dt/dx)*(
									MassFluxYFluxX(den_east, px_east, py_east,m_east,vel_snd)-
									MassFluxYFluxX(den_west, px_west, py_west,m_west,vel_snd))
									-(dt/dy)*(
									MassFluxYFluxY(den_north, px_north, py_north,m_north,vel_snd)-
									MassFluxYFluxY(den_south, px_south, py_south,m_south,vel_snd));
			}
		}	
	
}
		


void Fluid2D::CFLCondition(){ 
		dx = lengX / ( float ) ( Nx - 1 );
		dy = lengY / ( float ) ( Ny - 1 );
		dt = dx/10.0f;
}



float  Fluid2D::DensityFluxX(float n, float flx_x, float flx_y, float mass, float s){
	float f_1;
	f_1 = flx_x;
	return f_1;
}
float  Fluid2D::DensityFluxY(float n,float flx_x, float flx_y,float mass, float s){
	float f_1;
	f_1 = flx_y;
	return f_1;
}
float  Fluid2D::DensitySource(float n, float vel_x, float vel_y, float s){
	float q_1 =0.0f;
	return q_1;
}
float  Fluid2D::MassFluxXFluxX(float n,float flx_x, float flx_y,float mass, float s){
	float f_2;
	f_2 = flx_x * flx_x / n + n;
	return f_2;
}
float  Fluid2D::MassFluxXFluxY(float n,float flx_x, float flx_y,float mass, float s){
	float f_2;
	f_2 = flx_x * flx_y / n;
	return f_2;
}
float  Fluid2D::MassFluxYFluxX(float n,float flx_x, float flx_y,float mass, float s){
	float f_3;
	f_3 = flx_x * flx_y / n;
	return f_3;
}
float  Fluid2D::MassFluxYFluxY(float n,float flx_x, float flx_y,float mass, float s){
	float f_3;
	f_3 = flx_y * flx_y / n + n;
	return f_3;
}
float  Fluid2D::MassFluxXSource(float n,float flx_x, float flx_y, float s){
	float q_2 =0.0f;
	return q_2;
}
float  Fluid2D::MassFluxYSource(float n,float flx_x, float flx_y, float s){
	float q_3 =0.0f;
	return q_3;
}	

void Fluid2D::SetFileName(){
	char buffer [50];
	sprintf (buffer, "S=%.2fvis=%.2f", vel_snd, kin_vis);
	file_infix = buffer;
}

void Fluid2D::CreateFluidFile(){
	this->SetFileName();
	std::string previewfile = "preview_2D_" + file_infix + ".dat" ;
	data_preview.open (previewfile);
	data_preview << scientific; 
}

void Fluid2D::WriteFluidFile(float t){
int j=Ny/2;
data_preview << t << "\t" << Den[Nx - 1 + j * Nx] << "\t" << FlxX[Nx - 1 + j * Nx] << "\t" << Den[0 + j * Nx] << "\t" << FlxX[0 + j * Nx] << "\n";
}



void Fluid2D::SetSimulationTime(){
	Tmax=5.0f+0.02f*vel_snd+20.0f/vel_snd;
}

GrapheneFluid2D::GrapheneFluid2D(int size_nx, int size_ny, float sound_velocity, float fermi_velocity, float shear_viscosity, float collision_frequency): Fluid2D(size_nx, size_ny, sound_velocity,  shear_viscosity){
	vel_fer =fermi_velocity;
	col_freq =collision_frequency;
	char buffer [50];
	sprintf (buffer, "south=%.2fvF=%.2fvis=%.2fl=%.2f", vel_snd, vel_fer, kin_vis, col_freq);
	file_infix = buffer;
}


void GrapheneFluid2D::SetSimulationTime(){
	float s;
	s=this->GetVelSnd();
	this->SetTmax(5.0f+0.02f*s+20.0f/s);
}

void GrapheneFluid2D::MassFluxToVelocity(){
	for(int c=0; c <= Nx * Ny - 1; c++){
		VelX[c]= FlxX[c] * pow(Den[c], -1.5f);
		VelY[c]= FlxY[c] * pow(Den[c], -1.5f);
		CurX[c] = VelX[c] * Den[c];
		CurY[c] = VelY[c] * Den[c];
	}
}


void GrapheneFluid2D::SetVelFer(float x){ vel_fer=x;}
float GrapheneFluid2D::GetVelFer(){ return vel_fer;  }
void GrapheneFluid2D::SetColFreq(float x){ col_freq=x; }
float GrapheneFluid2D::GetColFreq(){ return col_freq; }


void GrapheneFluid2D::CFLCondition(){ // Eventual redefinition 
	dx = lengX / ( float ) ( Nx - 1 );
	dy = lengY / ( float ) ( Ny - 1 );
	//dt = 2.4/(vel_snd*sqrt(25.0/(dx*dx)+16.0/(dy*dy)));
	float lambda;
	if(vel_snd<0.36f*vel_fer){
		lambda=1.2f*vel_fer;
	}else{
		lambda=1.97f*vel_snd + 0.5f*vel_fer;
	}
	dt = dx/lambda;
}	

float  GrapheneFluid2D::DensityFluxX(float n,float flx_x, float flx_y,float mass, float s){ // Double Please Review this João
	float f_1;
	f_1 = flx_x / sqrt(n);
	return f_1;
}
float  GrapheneFluid2D::DensityFluxY(float n,float flx_x, float flx_y,float mass, float s){
	float f_1;
	f_1 = flx_y / sqrt(n);
	return f_1;
}
// 27% of cpu usage of them 22.6% are int the pow function
float  GrapheneFluid2D::MassFluxXFluxX(float n,float flx_x, float flx_y,float mass, float s){
	float f_2;
	f_2 = flx_x * flx_x / mass + vel_fer * vel_fer * mass / 3.0f + 0.5f * s * s * n * n;
	return f_2;
}
float  GrapheneFluid2D::MassFluxXFluxY(float n,float flx_x, float flx_y,float mass, float s){
	float f_2;
	f_2 = flx_x * flx_y / mass;
	return f_2;
}
float  GrapheneFluid2D::MassFluxYFluxX(float n,float flx_x, float flx_y,float mass, float s){
	float f_3;
	f_3 = flx_x * flx_y / mass;
	return f_3;
}
float  GrapheneFluid2D::MassFluxYFluxY(float n,float flx_x, float flx_y,float mass, float s){
	float f_3;
	f_3 = flx_y * flx_y / mass + vel_fer * vel_fer * mass / 3.0f + 0.5f * s * s * n * n;
	return f_3;
}




// Pedro: para ja nao vamos incluir sources 
float  GrapheneFluid2D::DensitySource(float n,float flx_x, float flx_y, float s){
	float q_1 =0;
	return q_1;
}
float  GrapheneFluid2D::MassFluxXSource(float n,float flx_x, float flx_y, float s){
	float q_2 =0;
	return q_2;
}
float  GrapheneFluid2D::MassFluxYSource(float n,float flx_x, float flx_y, float s){
	float q_3 =0;
	return q_3;
}	


void GrapheneFluid2D::MagneticSource(){
	float px_0,py_0,sqrtn_0;
	float wc=10.0;
	for(int kp=1+Nx; kp<=Nx*Ny-Nx-2; kp++){ //correr a grelha principal evitando as fronteiras
		if( kp%Nx!=Nx-1 && kp%Nx!=0){
			sqrtn_0=sqrt(Den[kp]);
			px_0=FlxX[kp];
			py_0=FlxY[kp];
			FlxX[kp]= px_0 * cos(wc * dt / sqrtn_0) - py_0 * sin(wc * dt / sqrtn_0);
			FlxY[kp]= px_0 * sin(wc * dt / sqrtn_0) + py_0 * cos(wc * dt / sqrtn_0);
		}
	}		
}

void GrapheneFluid2D::ViscosityFTCS() {
int north, south, east, west;
float mass_den_center, mass_den_north, mass_den_south, mass_den_east, mass_den_west;
//calculate laplacians
	for (int kp = 1 + Nx; kp <= Nx * Ny - Nx - 2; kp++) { //correr a grelha principal evitando as fronteiras
		div_t divresult;
		divresult = div(kp, Nx);
		int j = divresult.quot;
		int i = divresult.rem;
		if (kp % Nx != Nx - 1 && kp % Nx != 0) {
			north = i + (j + 1) * Nx;
			south = i + (j - 1) * Nx;
			east = i + 1 + j * Nx;
			west = i - 1 + j * Nx;
			mass_den_center = pow(Den[kp], 1.5f);
			mass_den_north = pow(Den[north], 1.5f);
			mass_den_south = pow(Den[south], 1.5f);
			mass_den_east = pow(Den[east], 1.5f);
			mass_den_west = pow(Den[west], 1.5f);
			lap_flxX[kp] =
					(-4.0f * FlxX[kp] / mass_den_center + FlxX[north] / mass_den_north + FlxX[south] / mass_den_south + FlxX[east] / mass_den_east +
					 FlxX[west] / mass_den_west) / (dx * dx);
			lap_flxY[kp] =
					(-4.0f * FlxY[kp] / mass_den_center + FlxY[north] / mass_den_north + FlxY[south] / mass_den_south + FlxY[east] / mass_den_east +
					 FlxY[west] / mass_den_west) / (dx * dx);
		}
	}

	//FTCS algorithm
	float old_px,old_py;
	for (int kp = 1 + Nx; kp <= Nx * Ny - Nx - 2; kp++) { //correr a grelha principal evitando as fronteiras
		if (kp % Nx != Nx - 1 && kp % Nx != 0) {
			old_px=FlxX[kp];
			old_py=FlxY[kp];
			FlxX[kp] = old_px + dt * (kin_vis * lap_flxX[kp] );
			FlxY[kp] = old_py + dt * (kin_vis * lap_flxY[kp] );
		}
	}
}

void GrapheneFluid2D::SourceFTCS(){
	float px_0,py_0,sqrtn_0;
	float wc=10.0;
	for(int kp=1+Nx; kp<=Nx*Ny-Nx-2; kp++){ //correr a grelha principal evitando as fronteiras
		if( kp%Nx!=Nx-1 && kp%Nx!=0){
			sqrtn_0=sqrt(Den[kp]);
			px_0=FlxX[kp];
			py_0=FlxY[kp];
			FlxX[kp]= px_0 - 0.5f * dt * wc * py_0 / sqrtn_0;
			FlxY[kp]= py_0 + 0.5f * dt * wc * px_0 / sqrtn_0;
		}
	}
}

void GrapheneFluid2D::WriteAtributes(){
	const FloatType      hdf5_float(PredType::NATIVE_FLOAT);
	const IntType        hdf5_int(PredType::NATIVE_INT);
	int total_steps= static_cast<int>(Tmax / dt);
	//Create the data space for the attribute.
	hsize_t dim_atr[1] = { 1 };
	DataSpace atr_dataspace = DataSpace (1, dim_atr );
	// Create a group attribute. 
	Attribute atr_vel_snd  = GrpDat->createAttribute("south parameter", hdf5_float, atr_dataspace);
	Attribute atr_vel_fer  = GrpDat->createAttribute("Fermi velocity", hdf5_float, atr_dataspace);
	Attribute atr_kin_vis = GrpDat->createAttribute("Kinetic viscosity", hdf5_float, atr_dataspace);
	Attribute atr_col_freq = GrpDat->createAttribute("Collision frequency", hdf5_float, atr_dataspace);
	Attribute atr_dx = GrpDat->createAttribute("Space discretisation step", hdf5_float, atr_dataspace);
	Attribute atr_dt = GrpDat->createAttribute("Time discretisation step", hdf5_float, atr_dataspace);
	Attribute atr_total_time = GrpDat->createAttribute("Total simulation time", hdf5_float, atr_dataspace);
	Attribute atr_num_space_points = GrpDat->createAttribute("Number of spatial points", hdf5_int, atr_dataspace);
	Attribute atr_num_time_steps = GrpDat->createAttribute("Number of time steps", hdf5_int, atr_dataspace);
	// Write the attribute data.
	atr_vel_snd.write( hdf5_float, &vel_snd);
	atr_vel_fer.write( hdf5_float, &vel_fer);
	atr_col_freq.write(hdf5_float, &col_freq);
	atr_kin_vis.write(hdf5_float, &kin_vis); 
	atr_dx.write(hdf5_float, &dx);
	atr_dt.write( hdf5_float, &dt);
	atr_num_space_points.write( hdf5_int, &Nx);
	atr_total_time.write( hdf5_float, &Tmax);
	atr_num_time_steps.write(hdf5_int, &total_steps);
	// Close the attributes.
	atr_num_time_steps.close();
	atr_col_freq.close();
	atr_vel_fer.close();
	atr_vel_snd.close();
	atr_kin_vis.close();
	atr_dx.close();
	atr_dt.close();
	atr_total_time.close();
	atr_num_space_points.close();
}
