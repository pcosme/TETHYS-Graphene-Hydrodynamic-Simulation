// 2D version



#include "includes/Fluid2DLib.h"
#include "SetUpParametersLib.h"
#include "GrapheneFluid2DLib.h"


using namespace H5;
using namespace std;


Fluid2D::Fluid2D(const SetUpParameters &input_parameters) : TethysBase{input_parameters.SizeX, input_parameters.SizeY, 2}{

	Nx = input_parameters.SizeX;
	Ny = input_parameters.SizeY;
	lengX=input_parameters.Length;
	lengY=input_parameters.Width;
	dx = lengX / ( float ) ( Nx - 1 );
	dy = lengY / ( float ) ( Ny - 1 );
	vel_snd = input_parameters.SoundVelocity;//sound_velocity;
	kin_vis = input_parameters.ShearViscosity;//shear_viscosity;

	char buffer [50];
	sprintf (buffer, "S=%.2fvF=%.2fvis=%.2fl=%.2fwc=%.2f", vel_snd, vel_fer, kin_vis, col_freq,cyc_freq);
	file_infix = buffer;
	// main grid variables Nx*Ny
	Den 		= new float[Nx * Ny]();
	VelX 		= new float[Nx * Ny]();
	VelY 		= new float[Nx * Ny]();
	FlxX 		= new float[Nx * Ny]();
	FlxY 		= new float[Nx * Ny]();
	CurX 		= new float[Nx * Ny]();
	CurY 		= new float[Nx * Ny]();
	vel_snd_arr	= new float[Nx * Ny]();
	velX_dx	= new float[Nx * Ny]();
	velX_dy	= new float[Nx * Ny]();
	velY_dx	= new float[Nx * Ny]();
	velY_dy	= new float[Nx * Ny]();

	lap_flxX = new float[Nx*Ny](); //new grids for the laplacians
	lap_flxY = new float[Nx*Ny](); //in fact they could be smaller but thiw way they are just 0 at the borders who do not evolve

	// 1st Aux. Grid variables (Nx-1)*(Ny-1)
	den_mid		= new float[(Nx-1)*(Ny-1)]();
	velX_dx_mid	= new float[(Nx-1)*(Ny-1)]();
	velX_dy_mid	= new float[(Nx-1)*(Ny-1)]();
	velY_dx_mid	= new float[(Nx-1)*(Ny-1)]();
	velY_dy_mid	= new float[(Nx-1)*(Ny-1)]();
	flxX_mid	= new float[(Nx-1)*(Ny-1)]();
	flxY_mid	= new float[(Nx-1)*(Ny-1)]();
	vel_snd_arr_mid	= new float[(Nx-1)*(Ny-1)]();

}

Fluid2D::~Fluid2D() = default;

void Fluid2D::SetSound(){
	for(int kp=0; kp<=Nx*Ny-1; kp++) { //correr a grelha principal evitando as fronteiras
		div_t divresult;
		divresult = div(kp, Nx);
		auto j = static_cast<float>(divresult.quot);
		auto i = static_cast<float>(divresult.rem);
		vel_snd_arr[kp]= Sound_Velocity_Anisotropy(i*dx, j*dy , vel_snd);
	}
	for(int ks=0; ks<=Nx*Ny-Nx-Ny; ks++) { //correr todos os pontos da grelha secundaria
		div_t divresult;
		divresult = div(ks, Nx - 1);
		auto j = static_cast<float>(divresult.quot);
		auto i = static_cast<float>(divresult.rem);
		vel_snd_arr_mid[ks]= Sound_Velocity_Anisotropy((i+0.5f)*dx, (j+0.5f)*dy , vel_snd);
	}
}




void Fluid2D::InitialCondRand(){
	random_device rd;
	float maxrand;
	maxrand = (float) random_device::max();
//#pragma omp parallel for default(none) shared(Nx,Ny,maxrand,rd)
	for (int c = 0; c < Nx*Ny; c++ ){
		float noise;
		noise =  (float) rd()/maxrand ; //(float) rand()/ (float) RAND_MAX ;
		Den[c] = 1.0f + 0.005f * (noise - 0.5f);
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
//#pragma omp parallel for default(none) shared(VelX,VelY,FlxX,FlxY,Den)
	for(int c=0; c <= Nx * Ny - 1; c++){
		VelX[c]= FlxX[c] / Den[c];
		VelY[c]= FlxY[c] / Den[c];
	}
}

void Fluid2D::VelocityToCurrent() {
//#pragma omp parallel for default(none) shared(Nx,Ny,VelX,VelY,CurX,CurY,Den)
	for(int c=0; c <= Nx * Ny - 1; c++){
		CurX[c] = VelX[c] * Den[c];
		CurY[c] = VelY[c] * Den[c];
	}
}

void Fluid2D::Richtmyer(){

	this->VelocityGradient();

#pragma omp parallel for default(none) shared(Nx,Ny,FlxX,FlxY,Den,flxX_mid,flxY_mid,den_mid,vel_snd_arr,dt,dx)
		for(int ks=0; ks<=Nx*Ny-Nx-Ny; ks++){ //correr todos os pontos da grelha secundaria de den_mid
			int northeast,northwest,southeast,southwest;
			float den_north, den_south ,den_east ,den_west, px_north, px_south, px_east, px_west, py_north, py_south, py_east, py_west,m_east,m_west,m_north,m_south;
			float  sound_north, sound_south ,sound_east ,sound_west;


			GridPoint point(ks,Nx,Ny,true);
			northeast=point.NE;
			northwest=point.NW;
			southeast=point.SE;
			southwest=point.SW;


			sound_north = 0.5f * (vel_snd_arr[northeast] + vel_snd_arr[northwest]);
			sound_south = 0.5f*(vel_snd_arr[southeast] + vel_snd_arr[southwest]);
			sound_east = 0.5f*(vel_snd_arr[northeast] + vel_snd_arr[southeast]);
			sound_west = 0.5f*(vel_snd_arr[northwest] + vel_snd_arr[southwest]);

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

			m_east=sqrt(den_east*den_east*den_east);
			m_west=sqrt(den_west*den_west*den_west);
			m_north=sqrt(den_north*den_north*den_north);
			m_south=sqrt(den_south*den_south*den_south);

			float den_avg = 0.25f * (Den[southwest] + Den[southeast] + Den[northwest] + Den[northeast]);
			float flx_x_avg = 0.25f * (FlxX[southwest] + FlxX[southeast] + FlxX[northwest] + FlxX[northeast]);
			float flx_y_avg = 0.25f * (FlxY[southwest] + FlxY[southeast] + FlxY[northwest] + FlxY[northeast]);

			den_mid[ks] = den_avg
					-0.5f*(dt/dx)*(
						DensityFluxX(den_east, px_east, py_east,m_east,sound_east)-
						DensityFluxX(den_west, px_west, py_west,m_west,sound_west))
					-0.5f*(dt/dy)*(
						DensityFluxY(den_north, px_north, py_north,m_north,sound_north)-
						DensityFluxY(den_south, px_south, py_south,m_south,sound_south))//;
					+0.5f*dt*DensitySource(den_avg, flx_x_avg, flx_y_avg, 0.0f, 0.0f);


			flxX_mid[ks] = flx_x_avg
					-0.5f*(dt/dx)*(
						MassFluxXFluxX(den_east, px_east, py_east,m_east,sound_east)-
						MassFluxXFluxX(den_west, px_west, py_west,m_west,sound_west))
					-0.5f*(dt/dy)*(
						MassFluxXFluxY(den_north, px_north, py_north,m_north,sound_north)-
						MassFluxXFluxY(den_south, px_south, py_south,m_south,sound_south))//;
					+0.5f*dt*MassFluxXSource(den_avg, flx_x_avg, flx_y_avg, 0.0f, 0.0f);

			flxY_mid[ks] = flx_y_avg
					-0.5f*(dt/dx)*(
						MassFluxYFluxX(den_east, px_east, py_east,m_east,sound_east)-
						MassFluxYFluxX(den_west, px_west, py_west,m_west,sound_west))
					-0.5f*(dt/dy)*(
						MassFluxYFluxY(den_north, px_north, py_north,m_north,sound_north)-
						MassFluxYFluxY(den_south, px_south, py_south,m_south,sound_south))//;
					+0.5f*dt*MassFluxYSource(den_avg, flx_x_avg, flx_y_avg, 0.0f, 0.0f);
		}

	this->VelocityGradientMid();

#pragma omp parallel for default(none) shared(Nx,Ny,FlxX,FlxY,Den,flxX_mid,flxY_mid,den_mid,vel_snd_arr_mid,dt,dx)
		for(int kp=1+Nx; kp<=Nx*Ny-Nx-2; kp++){ //correr a grelha principal evitando as fronteiras
			int northeast,northwest,southeast,southwest;
			float den_north, den_south ,den_east ,den_west, px_north, px_south, px_east, px_west, py_north, py_south, py_east, py_west,m_east,m_west,m_north,m_south;
			float  sound_north, sound_south ,sound_east ,sound_west;

			GridPoint point(kp,Nx,Ny,false);

			if( kp%Nx!=Nx-1 && kp%Nx!=0){

				northeast=point.NE;
				northwest=point.NW;
				southeast=point.SE;
				southwest=point.SW;

				sound_north = 0.5f * (vel_snd_arr_mid[northeast] + vel_snd_arr_mid[northwest]);
				sound_south = 0.5f*(vel_snd_arr_mid[southeast] + vel_snd_arr_mid[southwest]);
				sound_east = 0.5f*(vel_snd_arr_mid[northeast] + vel_snd_arr_mid[southeast]);
				sound_west = 0.5f*(vel_snd_arr_mid[northwest] + vel_snd_arr_mid[southwest]);

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

				m_east=sqrt(den_east*den_east*den_east);
				m_west=sqrt(den_west*den_west*den_west);
				m_north=sqrt(den_north*den_north*den_north);
				m_south=sqrt(den_south*den_south*den_south);

				float den_old = Den[kp];
				float flx_x_old = FlxX[kp];
				float flx_y_old = FlxY[kp];


				Den[kp] = den_old - (dt / dx) * (
							DensityFluxX(den_east, px_east, py_east,m_east,sound_east)-
							DensityFluxX(den_west, px_west, py_west,m_west,sound_west))
						-(dt/dy)*(
							DensityFluxY(den_north, px_north, py_north,m_north,sound_north)-
							DensityFluxY(den_south, px_south, py_south,m_south,sound_south))
						+dt*DensitySource(den_old, flx_x_old, flx_y_old, 0.0f, 0.0f);
				FlxX[kp] = flx_x_old - (dt / dx) * (
							MassFluxXFluxX(den_east, px_east, py_east,m_east,sound_east)-
							MassFluxXFluxX(den_west, px_west, py_west,m_west,sound_west))
						-(dt/dy)*(
							MassFluxXFluxY(den_north, px_north, py_north,m_north,sound_north)-
							MassFluxXFluxY(den_south, px_south, py_south,m_south,sound_south))
						+dt*MassFluxXSource(den_old, flx_x_old, flx_y_old, 0.0f, 0.0f);
				FlxY[kp] = flx_y_old - (dt / dx) * (
							MassFluxYFluxX(den_east, px_east, py_east,m_east,sound_east)-
							MassFluxYFluxX(den_west, px_west, py_west,m_west,sound_west))
						-(dt/dy)*(
							MassFluxYFluxY(den_north, px_north, py_north,m_north,sound_north)-
							MassFluxYFluxY(den_south, px_south, py_south,m_south,sound_south))
						+dt*MassFluxYSource(den_old, flx_x_old, flx_y_old, 0.0f, 0.0f);

			}
		}
}

void Fluid2D::CflCondition(){
		dx = lengX / ( float ) ( Nx - 1 );
		dy = lengY / ( float ) ( Ny - 1 );
		dt = dx/10.0f;
}



float  Fluid2D::DensityFluxX(__attribute__((unused)) float n, float flx_x, __attribute__((unused)) float flx_y, __attribute__((unused)) float mass, __attribute__((unused)) float s){
	float f_1;
	f_1 = flx_x;
	return f_1;
}
float  Fluid2D::DensityFluxY(__attribute__((unused)) float n, __attribute__((unused)) float flx_x, float flx_y, __attribute__((unused)) float mass, __attribute__((unused)) float s){
	float f_1;
	f_1 = flx_y;
	return f_1;
}


float  Fluid2D::MassFluxXFluxX(float n,float flx_x,__attribute__((unused)) float flx_y,__attribute__((unused)) float mass,__attribute__((unused))  float s){
	float f_2;
	f_2 = flx_x * flx_x / n + n;
	return f_2;
}
float  Fluid2D::MassFluxXFluxY(float n,float flx_x, float flx_y,__attribute__((unused)) float mass,__attribute__((unused))  float s){
	float f_2;
	f_2 = flx_x * flx_y / n;
	return f_2;
}
float  Fluid2D::MassFluxYFluxX(float n,float flx_x, float flx_y,__attribute__((unused)) float mass,__attribute__((unused))  float s){
	float f_3;
	f_3 = flx_x * flx_y / n;
	return f_3;
}
float  Fluid2D::MassFluxYFluxY(float n, __attribute__((unused)) float flx_x, float flx_y,__attribute__((unused)) float mass,__attribute__((unused))  float s){
	float f_3;
	f_3 = flx_y * flx_y / n + n;
	return f_3;
}


void Fluid2D::CreateFluidFile(){
	std::string previewfile = "preview_2D_" + file_infix + ".dat" ;
	data_preview.open (previewfile);
	data_preview << scientific; 
}

void Fluid2D::WriteFluidFile(float t){
	int j=Ny/2;
	int pos_end = Nx - 1 + j*Nx ;
	int pos_ini = j*Nx ;
		if(!isfinite(Den[pos_end]) || !isfinite(Den[pos_ini]) || !isfinite(FlxX[pos_end]) || !isfinite(FlxX[pos_ini])){
			cerr << "ERROR: numerical method failed to converge" <<"\nExiting"<< endl;
			CloseHdf5File();
			exit(EXIT_FAILURE);
		}
	data_preview << t << "\t"
	<< Den[pos_end]  << "\t"
	<< FlxX[pos_end] << "\t"
	<< Den[pos_ini]  << "\t"
	<< FlxX[pos_ini] << "\n";
}

void Fluid2D::SetSimulationTime(){
	Tmax=5.0f+0.02f*vel_snd+20.0f/vel_snd;
}

/*
void Fluid2D::VelocityLaplacian() {
	this->MassFluxToVelocity();

//calculate laplacians
//#pragma omp parallel for  default(none) shared(lap_flxX,lap_flxY,VelX,VelY,Nx,Ny)

	for(int ks=0; ks<=Nx*Ny-Nx-Ny; ks++) { //correr todos os pontos da grelha secundaria de den_mid
		int northeast, northwest, southeast, southwest;
		div_t divresult;
		divresult = div(ks, Nx - 1);
		int j = divresult.quot;
		int i = divresult.rem;
		northeast = i + 1 + (j + 1) * Nx;
		northwest = i + (j + 1) * Nx;
		southeast = i + 1 + j * Nx;
		southwest = i + j * Nx;

		velX_lap_mid[ks] = (-4.0f * velX_mid[ks]  + VelX[northeast]  + VelX[southeast]  + VelX[northwest]  +
		                    VelX[southwest] )/(0.5f*dx*dx);
		velY_lap_mid[ks] = (-4.0f * velY_mid[ks]  + VelY[northeast]  + VelY[southeast]  + VelY[northwest]  +
		                    VelY[southwest] )/(0.5f*dx*dx);
	}


	for(int kp=1+Nx; kp<=Nx*Ny-Nx-2; kp++){ //correr a grelha principal evitando as fronteiras
		int northeast,northwest,southeast,southwest;

		div_t divresult;
		divresult = div (kp,Nx);
		int j=divresult.quot;
		int i=divresult.rem;

		if( kp%Nx!=Nx-1 && kp%Nx!=0) {
			northeast = i + j * (Nx - 1);
			northwest = i - 1 + j * (Nx - 1);
			southeast = i + (j - 1) * (Nx - 1);
			southwest = i - 1 + (j - 1) * (Nx - 1);

			velX_lap[kp] = (-4.0f * VelX[kp]  + velX_mid[northeast]  + velX_mid[southeast]  + velX_mid[northwest]  +
					velX_mid[southwest] )/(0.5f*dx*dx);
			velY_lap[kp] = (-4.0f * VelY[kp]  + velY_mid[northeast]  + velY_mid[southeast]  + velY_mid[northwest]  +
			                velY_mid[southwest] )/(0.5f*dx*dx);

		}
	}
}
*/


void Fluid2D::VelocityLaplacianFtcs() {
	this->MassFluxToVelocity();

//calculate laplacians
#pragma omp parallel for  default(none) shared(lap_flxX,lap_flxY,VelX,VelY,Nx,Ny)
	for (int kp = 1 + Nx; kp <= Nx * Ny - Nx - 2; kp++) { //correr a grelha principal evitando as fronteiras
		int north, south, east, west;
		div_t divresult;
		divresult = div(kp, Nx);
		int j;
		j = divresult.quot;
		int i;
		i = divresult.rem;
		if (kp % Nx != Nx - 1 && kp % Nx != 0){
			north = i + (j + 1) * Nx;
			south = i + (j - 1) * Nx;
			east = i + 1 + j * Nx;
			west = i - 1 + j * Nx;
			lap_flxX[kp] =
					kin_vis*dt*(-4.0f * VelX[kp]  + VelX[north]  + VelX[south]  + VelX[east]  +
							VelX[west] ) / (dx * dx);
			lap_flxY[kp] =
					kin_vis*dt*(-4.0f * VelY[kp]  + VelY[north]  + VelY[south]  + VelY[east]  +
							VelY[west] ) / (dx * dx);
		}
	}
}

void Fluid2D::VelocityLaplacianWeighted19() {
	this->MassFluxToVelocity();
	float sx=kin_vis*dt/(dx*dx);
	float sy=kin_vis*dt/(dy*dy);
#pragma omp parallel for default(none) shared(lap_flxX,lap_flxY,VelX,VelY,sx,sy)
	for (int kp = 1 + Nx; kp <= Nx * Ny - Nx - 2; kp++) { //correr a grelha principal evitando as fronteiras
		int north, south, east, west, northeast, northwest,southeast, southwest ;
		div_t divresult;
		divresult = div(kp, Nx);
		int j;
		j = divresult.quot;
		int i;
		i = divresult.rem;
		if (kp % Nx != Nx - 1 && kp % Nx != 0){
			north = i + (j + 1) * Nx;
			south = i + (j - 1) * Nx;
			east = i + 1 + j * Nx;
			west = i - 1 + j * Nx;
			northeast = i + 1 + (j + 1) * Nx;
			northwest = i - 1 + (j + 1) * Nx;
			southeast = i + 1 + (j - 1) * Nx;
			southwest = i - 1 + (j - 1) * Nx;
			lap_flxX[kp] = (4.0f*sx*sy-2.0f*sx-2.0f*sy)*VelX[kp] +
			               sx*sy*( VelX[northeast] + VelX[southeast] + VelX[northwest] + VelX[southwest])
			               + sy*(1.0f-2.0f*sx)*(VelX[north] + VelX[south])
			               + sx*(1.0f-2.0f*sy)*(VelX[west] + VelX[east]);
			lap_flxY[kp] = (4.0f*sx*sy-2.0f*sx-2.0f*sy)*VelY[kp] +
			               sx*sy*( VelY[northeast] + VelY[southeast] + VelY[northwest] + VelY[southwest])
			               + sy*(1.0f-2.0f*sx)*(VelY[north] + VelY[south])
			               + sx*(1.0f-2.0f*sy)*(VelY[west] + VelY[east]);
		}
	}
}


void Fluid2D:: ParabolicOperatorFtcs() {
	VelocityLaplacianFtcs();
	ForwardTimeOperator();
}

void Fluid2D::ParabolicOperatorWeightedExplicit19() {
	VelocityLaplacianWeighted19();
	ForwardTimeOperator();
}


void Fluid2D::SaveSound() {
	DataSet dataset_vel_snd = GrpDat->createDataSet("Sound velocicity", HDF5FLOAT, *DataspaceVelSnd);
	dataset_vel_snd.write(vel_snd_arr, HDF5FLOAT);
	dataset_vel_snd.close();
}

void Fluid2D::ReadSnapShot(const H5std_string &snap_name) {
	DataSet dataset_den = GrpDen->openDataSet( snap_name );
	DataSet dataset_vel_x  = GrpVelX->openDataSet( snap_name );
	DataSet dataset_vel_y = GrpVelY->openDataSet( snap_name );
	DataSpace dataspace = dataset_den.getSpace(); // Get dataspace of the dataset.
	if(dataset_den.attrExists ("time")){
		Attribute attr_time = dataset_den.openAttribute("time");
		attr_time.read(attr_time.getDataType(), &TimeStamp);
	}else{
		TimeStamp+=1.0f;
	}
	if(dataset_den.attrExists ("time step")){
		Attribute attr_counter = dataset_den.openAttribute("time step");
		attr_counter.read(attr_counter.getDataType(), &TimeStepCounter);
	}else{
		TimeStepCounter++;
	}
	dataset_den.read( Den, PredType::NATIVE_FLOAT,   dataspace );
	dataset_vel_x.read( VelX, PredType::NATIVE_FLOAT,   dataspace );
	dataset_vel_y.read( VelY, PredType::NATIVE_FLOAT,   dataspace );
	dataset_den.close();
	dataset_vel_x.close();
	dataset_vel_y.close();
}

void Fluid2D::SaveSnapShot() {
	hsize_t dim_atr[1] = { 1 };
	DataSpace atr_dataspace = DataSpace (1, dim_atr );

	int points_per_period = static_cast<int>((2.0 * MAT_PI / this->RealFreq()) / dt);
	snapshot_step = points_per_period / snapshot_per_period;

	this->MassFluxToVelocity();
	this->VelocityToCurrent();
	string str_time = to_string(TimeStepCounter / snapshot_step);
	str_time.insert(str_time.begin(), 5 - str_time.length(), '0');
	string name_dataset = "snapshot_" + str_time;

	DataSet dataset_den = GrpDen->createDataSet(name_dataset, HDF5FLOAT, *DataspaceDen);
	Attribute atr_step_den = dataset_den.createAttribute("time step", HDF5INT, atr_dataspace);
	Attribute atr_time_den = dataset_den.createAttribute("time", HDF5FLOAT, atr_dataspace);
	float currenttime=static_cast<float>(TimeStepCounter) * dt;
	atr_step_den.write(HDF5INT, &TimeStepCounter);
	atr_time_den.write(HDF5FLOAT , &currenttime);
	atr_step_den.close();
	atr_time_den.close();

	dataset_den.write(Den, HDF5FLOAT);
	dataset_den.close();

	DataSet dataset_vel_x = GrpVelX->createDataSet(name_dataset, HDF5FLOAT, *DataspaceVelX);
	Attribute atr_step_vel_x = dataset_vel_x.createAttribute("time step", HDF5INT, atr_dataspace);
	Attribute atr_time_vel_x = dataset_vel_x.createAttribute("time", HDF5FLOAT, atr_dataspace);
	dataset_vel_x.write(VelX, HDF5FLOAT);
	dataset_vel_x.close();
	atr_step_vel_x.write(HDF5INT, &TimeStepCounter);
	atr_time_vel_x.write(HDF5FLOAT , &currenttime);
	atr_step_vel_x.close();
	atr_time_vel_x.close();

	DataSet dataset_vel_y = GrpVelY->createDataSet(name_dataset, HDF5FLOAT, *DataspaceVelY);
	Attribute atr_step_vel_y = dataset_vel_y.createAttribute("time step", HDF5INT, atr_dataspace);
	Attribute atr_time_vel_y = dataset_vel_y.createAttribute("time", HDF5FLOAT, atr_dataspace);
	dataset_vel_y.write(VelY, HDF5FLOAT);
	dataset_vel_y.close();
	atr_step_vel_y.write(HDF5INT, &TimeStepCounter);
	atr_time_vel_y.write(HDF5FLOAT , &currenttime);
	atr_step_vel_y.close();
	atr_time_vel_y.close();
}

int Fluid2D::GetSnapshotStep() const { return snapshot_step;}
int Fluid2D::GetSnapshotFreq() const {return snapshot_per_period;}

bool Fluid2D::Snapshot() const {
	bool state;
	if(TimeStepCounter % snapshot_step == 0){
		state = true;
	}else{
		state = false;
	}
	return state;
}

float Fluid2D::DensitySource(__attribute__((unused)) float n, __attribute__((unused)) float flx_x,__attribute__((unused)) float flx_y,__attribute__((unused)) float mass,__attribute__((unused)) float s) {
	return 0.0f;
}

float Fluid2D::MassFluxXSource(__attribute__((unused)) float n,__attribute__((unused)) float flx_x,__attribute__((unused)) float flx_y,__attribute__((unused)) float mass,__attribute__((unused)) float s) {
	return 0.0f;
}
float Fluid2D::MassFluxYSource(__attribute__((unused)) float n,__attribute__((unused)) float flx_x,__attribute__((unused)) float flx_y,__attribute__((unused)) float mass,__attribute__((unused)) float s) {
	return 0.0f;
}

void Fluid2D::ForwardTimeOperator() {
#pragma omp parallel for default(none) shared(Nx,Ny,FlxX,FlxY,lap_flxX,lap_flxY,dt,cyc_freq)
	for (int kp = 1 + Nx; kp <= Nx * Ny - Nx - 2; kp++) { //correr a grelha principal evitando as fronteiras
		float flx_x_old,flx_y_old;
		if (kp % Nx != Nx - 1 && kp % Nx != 0) {
			flx_x_old=FlxX[kp];
			flx_y_old=FlxY[kp];
			FlxX[kp] = flx_x_old + lap_flxX[kp];
			FlxY[kp] = flx_y_old + lap_flxY[kp];
		}
	}
}

void Fluid2D::VelocityGradient() {

	int stride = Nx;

	float m_east,m_west,m_north,m_south;
	for(int kp=1+Nx; kp<=Nx*Ny-Nx-2; kp++){  //correr a grelha principal evitando as fronteiras
		if( kp%stride!=stride-1 && kp%stride!=0){ //aqui podem ser diferencas centradas quer em x quer em y
			GridPoint point(kp,Nx,Ny,false);
			m_east =  DensityToMass(Den[point.E]);
			m_west =  DensityToMass(Den[point.W]);
			m_north = DensityToMass(Den[point.N]);
			m_south = DensityToMass(Den[point.S]);
			velX_dx[kp] = ( (FlxX[point.E]/m_east) -(FlxX[point.W]/m_west) )/(2.0f*dx);
			velX_dy[kp] = ( (FlxX[point.N]/m_north)-(FlxX[point.S]/m_south) )/(2.0f*dy);
			velY_dx[kp] = ( (FlxY[point.E]/m_east) -(FlxY[point.W]/m_west) )/(2.0f*dx);
			velY_dy[kp] = ( (FlxY[point.N]/m_north)-(FlxY[point.S]/m_south) )/(2.0f*dy);
		}
	}
	for(int i=1 ; i<=Nx-2; i++){ // topo rede principal, ou seja j=(Ny - 1)
		float m_top,m_southsouth;
		int top= i + (Ny - 1) * stride;
		GridPoint point(top,Nx,Ny,false);
		int southsouth= i + (Ny - 3) * stride;
		m_top =   DensityToMass(Den[point.C]);
		m_east =  DensityToMass(Den[point.E]);
		m_west =  DensityToMass(Den[point.W]);
		m_south = DensityToMass(Den[point.S]);
		m_southsouth = DensityToMass(Den[southsouth]);
		velX_dx[top] = ( (FlxX[point.E]/m_east)-(FlxX[point.W]/m_west) )/(2.0f*dx); //OK
		velX_dy[top] = ( 3.0f*(FlxX[point.C]/m_top) -4.0f*(FlxX[point.S]/m_south) +1.0f*(FlxX[southsouth]/m_southsouth) )/(2.0f*dy); //backward finite difference
		velY_dx[top] = ( (FlxY[point.E]/m_east)-(FlxY[point.W]/m_west) )/(2.0f*dx); //OK
		velY_dy[top] = ( 3.0f*(FlxY[point.C]/m_top) -4.0f*(FlxY[point.S]/m_south) +1.0f*(FlxY[southsouth]/m_southsouth)  )/(2.0f*dy); //backward finite difference
	}
	for(int i=1 ; i<=Nx-2; i++){ // fundo rede principal, ou seja j=0
		float m_bottom,m_northnorth;
		int bottom=i; //i+0*nx
		GridPoint point(bottom,Nx,Ny,false);
		int northnorth=i+2*stride;
		m_bottom = DensityToMass(Den[point.C]);
		m_east =   DensityToMass(Den[point.E]);
		m_west =   DensityToMass(Den[point.W]);
		m_north =  DensityToMass(Den[point.N]);
		m_northnorth = DensityToMass(Den[northnorth]);
		velX_dx[bottom] = ( (FlxX[point.E]/m_east)-(FlxX[point.W]/m_west) )/(2.0f*dx);  //OK
		velX_dy[bottom] = ( -3.0f*(FlxX[point.C]/m_bottom) +4.0f*(FlxX[point.N]/m_north)-1.0f*(FlxX[northnorth]/m_northnorth) )/(2.0f*dy); //forward finite difference
		velY_dx[bottom] = ( (FlxY[point.E]/m_east)-(FlxY[point.W]/m_west) )/(2.0f*dx); //OK
		velY_dy[bottom] = ( -3.0f*(FlxY[point.C]/m_bottom) +4.0f*(FlxY[point.N]/m_north)-1.0f*(FlxY[northnorth]/m_northnorth) )/(2.0f*dy); //forward finite difference
	}
	for(int j=1; j<=Ny-2;j++){ //lado esquerdo da rede principal ou seja i=0
		float m_left, m_easteast;
		int left = 0 + j*stride;
		GridPoint point(left,Nx,Ny,false);
		int easteast = left + 2;
		m_left =  DensityToMass(Den[point.C]);
		m_north = DensityToMass(Den[point.N]);
		m_south = DensityToMass(Den[point.S]);
		m_east =  DensityToMass(Den[point.E]);
		m_easteast = DensityToMass(Den[easteast]);
		velX_dx[left] = ( -3.0f*(FlxX[point.C]/m_left) +4.0f*(FlxX[point.E]/m_east)-1.0f*(FlxX[easteast]/m_easteast) )/(2.0f*dx); //forward difference
		velX_dy[left] = ( (FlxX[point.N]/m_north)-(FlxX[point.S]/m_south) )/(2.0f*dy); //OK
		velY_dx[left] = ( -3.0f*(FlxY[point.C]/m_left) +4.0f*(FlxY[point.E]/m_east)-1.0f*(FlxY[easteast]/m_easteast) )/(2.0f*dx); //forward difference
		velY_dy[left] = ( (FlxY[point.N]/m_north)-(FlxY[point.S]/m_south) )/(2.0f*dy); //OK
	}
	for(int j=1; j<=Ny-2;j++){ //lado direito da rede principal ou seja i=(Nx-1)
		float m_rigth, m_westwest;
		int right = (Nx-1) + j*stride;
		GridPoint point(right,Nx,Ny,false);
		int westwest = right-2;
		m_rigth = DensityToMass(Den[point.C]);
		m_north = DensityToMass(Den[point.N]);
		m_south = DensityToMass(Den[point.S]);
		m_west =  DensityToMass(Den[point.W]);
		m_westwest = DensityToMass(Den[westwest]);
		velX_dx[right] = ( 3.0f*(FlxX[point.C]/m_rigth) -4.0f*(FlxX[point.W]/m_west) +1.0f*(FlxX[westwest]/m_westwest) )/(2.0f*dx); //backwar difference
		velX_dy[right] = ( (FlxX[point.N]/m_north)-(FlxX[point.S]/m_south) )/(2.0f*dy); //OK
		velY_dx[right] = ( 3.0f*(FlxY[point.C]/m_rigth) -4.0f*(FlxY[point.W]/m_west) +1.0f*(FlxY[westwest]/m_westwest) )/(2.0f*dx);
		velY_dy[right] = ( (FlxY[point.N]/m_north)-(FlxY[point.S]/m_south) )/(2.0f*dy); //OK
	}
	//os 4 cantos em que ambas a derivadas nao podem ser centradas
	int kp;
	// i=0 j=0 forward x forward y
	kp = 0 + 0*Nx;
	float m_0 =  DensityToMass(Den[kp]);
	float m_x1 = DensityToMass(Den[kp+1]);
	float m_x2 = DensityToMass(Den[kp+2]);
	float m_y1 = DensityToMass(Den[kp+Nx]);
	float m_y2 = DensityToMass(Den[kp+2*Nx]);
	velX_dx[kp] = ( -3.0f*(FlxX[kp]/m_0 ) +4.0f*(FlxX[kp+1]/m_x1 )-1.0f*(FlxX[kp+2]/m_x2 )  )/(2.0f*dx);
	velX_dy[kp] = ( -3.0f*(FlxX[kp]/m_0 ) +4.0f*(FlxX[kp+1*Nx]/m_y1 )-1.0f*(FlxX[kp+2*Nx]/m_y2 )  )/(2.0f*dy);
	velY_dx[kp] = ( -3.0f*(FlxY[kp]/m_0 ) +4.0f*(FlxY[kp+1]/m_x1 )-1.0f*(FlxY[kp+2]/m_x2 )  )/(2.0f*dx);
	velY_dy[kp] = ( -3.0f*(FlxY[kp]/m_0 ) +4.0f*(FlxY[kp+1*Nx]/m_y1 )-1.0f*(FlxY[kp+2*Nx]/m_y2 )  )/(2.0f*dy);

	// i=(Nx-1) j=0 backward x forward y
	kp = (Nx-1) + 0*Nx;
	m_0 =  DensityToMass(Den[kp]);
	m_x1 = DensityToMass(Den[kp-1]);
	m_x2 = DensityToMass(Den[kp-2]);
	m_y1 = DensityToMass(Den[kp+Nx]);
	m_y2 = DensityToMass(Den[kp+2*Nx]);
	velX_dx[kp] = (  3.0f*(FlxX[kp]/m_0 ) -4.0f*(FlxX[kp-1]/m_x1 )+1.0f*(FlxX[kp-2]/m_x2 )  )/(2.0f*dx);
	velX_dy[kp] = ( -3.0f*(FlxX[kp]/m_0 ) +4.0f*(FlxX[kp+Nx]/m_y1 )-1.0f*(FlxX[kp+2*Nx]/m_y2 )  )/(2.0f*dy);
	velY_dx[kp] = (  3.0f*(FlxY[kp]/m_0 ) -4.0f*(FlxY[kp-1]/m_x1 )+1.0f*(FlxY[kp-2]/m_x2 )  )/(2.0f*dx);
	velY_dy[kp] = ( -3.0f*(FlxY[kp]/m_0 ) +4.0f*(FlxY[kp+Nx]/m_y1 )-1.0f*(FlxY[kp+2*Nx]/m_y2 )  )/(2.0f*dy);

	// i=0 j=(Ny-1) forward x backward y
	kp = 0 + (Ny-1)*Nx;
	m_0 =  DensityToMass(Den[kp]);
	m_x1 = DensityToMass(Den[kp+1]);
	m_x2 = DensityToMass(Den[kp+2]);
	m_y1 = DensityToMass(Den[kp-Nx]);
	m_y2 = DensityToMass(Den[kp-2*Nx]);
	velX_dx[kp] = ( -3.0f*(FlxX[kp]/m_0 ) +4.0f*(FlxX[kp+1]/m_x1 )-1.0f*(FlxX[kp+2]/m_x1 )  )/(2.0f*dx);
	velX_dy[kp] = (  3.0f*(FlxX[kp]/m_0 ) -4.0f*(FlxX[kp-Nx]/m_y1 )+1.0f*(FlxX[kp-2*Nx]/m_y2 )  )/(2.0f*dy);
	velY_dx[kp] = ( -3.0f*(FlxY[kp]/m_0 ) +4.0f*(FlxY[kp+1]/m_x1 )-1.0f*(FlxY[kp+2]/m_x2 )  )/(2.0f*dx);
	velY_dy[kp] = (  3.0f*(FlxY[kp]/m_0 ) -4.0f*(FlxY[kp-Nx]/m_y1 )+1.0f*(FlxY[kp-2*Nx]/m_y2 )  )/(2.0f*dy);

	// i=(Nx-1) j=(Ny-1) backward x backward y
	kp = (Nx-1) + (Ny-1)*Nx;
	m_0 =  DensityToMass(Den[kp]);
	m_x1 = DensityToMass(Den[kp-1]);
	m_x2 = DensityToMass(Den[kp-2]);
	m_y1 = DensityToMass(Den[kp-Nx]);
	m_y2 = DensityToMass(Den[kp-2*Nx]);
	velX_dx[kp] = ( 3.0f*(FlxX[kp]/m_0 ) -4.0f*(FlxX[kp-1]/m_x1 )+1.0f*(FlxX[kp-2]/m_x2 ) )/(2.0f*dx);
	velX_dy[kp] = ( 3.0f*(FlxX[kp]/m_0 ) -4.0f*(FlxX[kp-Nx]/m_y1 )+1.0f*(FlxX[kp-2*Nx]/m_y2 ) )/(2.0f*dy);
	velY_dx[kp] = ( 3.0f*(FlxY[kp]/m_0 ) -4.0f*(FlxY[kp-1]/m_x1 )+1.0f*(FlxY[kp-2]/m_x2 ) )/(2.0f*dx);
	velY_dy[kp] = ( 3.0f*(FlxY[kp]/m_0 ) -4.0f*(FlxY[kp-Nx]/m_y1 )+1.0f*(FlxY[kp-2*Nx]/m_y2 ) )/(2.0f*dy);

}



void Fluid2D::VelocityGradientMid() {
	int north, south, east, west;
	float m_east,m_west,m_north,m_south;
	for(int ks=1+(Nx-1); ks<=(Nx-2)+(Ny-2)*(Nx-1); ks++){ //correr todos os pontos da grelha secundaria de _mid EVITANDO FRONTEIRAS
		if( ks%(Nx-1)!=Nx-2 && ks%(Nx-1)!=0) {
			north = ks - (Nx - 1);
			south = ks + (Nx - 1);
			east = ks + 1;
			west = ks - 1;
			m_east = DensityToMass(den_mid[east]);
			m_west = DensityToMass(den_mid[west]);
			m_north = DensityToMass(den_mid[north]);
			m_south = DensityToMass(den_mid[south]);
			velX_dx_mid[ks] = ( (flxX_mid[east]/m_east) - (flxX_mid[west]/m_west) ) / (2.0f * dx);
			velX_dy_mid[ks] = ( (flxX_mid[south]/m_south) - (flxX_mid[north]/m_north) ) / (2.0f * dy);
			velY_dx_mid[ks] = ( (flxY_mid[east]/m_east) - (flxY_mid[west]/m_west) ) / (2.0f * dx);
			velY_dy_mid[ks] = ( (flxY_mid[south]/m_south) - (flxY_mid[north]/m_north) ) / (2.0f * dy);
		}
	}
	for(int i=1 ; i<=(Nx-1)-2; i++){ // topo rede principal, ou seja j=((Ny-1) - 1)
		float m_top,m_southsouth;
		int top= i + ((Ny-1) - 1) * (Nx-1);
		east=top+1;
		west=top-1;
		south= i + ((Ny-1) - 2) * (Nx-1);
		int southsouth= i + ((Ny-1) - 3) * (Nx-1);
		m_top = DensityToMass(den_mid[top]);
		m_east = DensityToMass(den_mid[east]);
		m_west = DensityToMass(den_mid[west]);
		m_south = DensityToMass(den_mid[south]);
		m_southsouth = DensityToMass(den_mid[southsouth]);
		velX_dx_mid[top] = ( (flxX_mid[east]/m_east)-(flxX_mid[west]/m_west) )/(2.0f*dx); //OK
		velX_dy_mid[top] = ( 3.0f*(flxX_mid[top]/m_top) -4.0f*(flxX_mid[south]/m_south) +1.0f*(flxX_mid[southsouth]/m_southsouth) )/(2.0f*dy); //backward finite difference
		velY_dx_mid[top] = ( (flxY_mid[east]/m_east)-(flxY_mid[west]/m_west) )/(2.0f*dx); //OK
		velY_dy_mid[top] = ( 3.0f*(flxY_mid[top]/m_top) -4.0f*(flxY_mid[south]/m_south) +1.0f*(flxY_mid[southsouth]/m_southsouth)  )/(2.0f*dy); //backward finite difference
	}
	for(int i=1 ; i<=(Nx-1)-2; i++){ // fundo rede principal, ou seja j=0
		float m_bottom,m_northnorth;
		int bottom=i; //i+0*nx
		east=bottom+1;
		west=bottom-1;
		north=i+(Nx-1);
		int northnorth=i+2*(Nx-1);
		m_bottom = DensityToMass(den_mid[bottom]);
		m_east = DensityToMass(den_mid[east]);
		m_west = DensityToMass(den_mid[west]);
		m_north = DensityToMass(den_mid[north]);
		m_northnorth = DensityToMass(den_mid[northnorth]);
		velX_dx_mid[bottom] = ( (flxX_mid[east]/m_east)-(flxX_mid[west]/m_west) )/(2.0f*dx);  //OK
		velX_dy_mid[bottom] = ( -3.0f*(flxX_mid[bottom]/m_bottom) +4.0f*(flxX_mid[north]/m_north)-1.0f*(flxX_mid[northnorth]/m_northnorth) )/(2.0f*dy); //forward finite difference
		velY_dx_mid[bottom] = ( (flxY_mid[east]/m_east)-(flxY_mid[west]/m_west) )/(2.0f*dx); //OK
		velY_dy_mid[bottom] = ( -3.0f*(flxY_mid[bottom]/m_bottom) +4.0f*(flxY_mid[north]/m_north)-1.0f*(flxY_mid[northnorth]/m_northnorth) )/(2.0f*dy); //forward finite difference
	}
	for(int j=1; j<=(Ny-1)-2;j++){ //lado esquerdo da rede principal ou seja i=0
		float m_left, m_easteast;
		int left = 0 + j*(Nx-1);
		int easteast = left + 2;
		east =left+1;
		north=left+(Nx-1);
		south=left-(Nx-1);
		m_north = DensityToMass(den_mid[north]);
		m_south = DensityToMass(den_mid[south]);
		m_east = DensityToMass(den_mid[east]);
		m_easteast = DensityToMass(den_mid[easteast]);
		m_left= DensityToMass(den_mid[left]);
		velX_dx_mid[left] = ( -3.0f*(flxX_mid[left]/m_left) +4.0f*(flxX_mid[east]/m_east)-1.0f*(flxX_mid[easteast]/m_easteast)  )/(2.0f*dx); //forward difference
		velX_dy_mid[left] = ( (flxX_mid[north]/m_north)-(flxX_mid[south]/m_south) )/(2.0f*dy); //OK
		velY_dx_mid[left] = ( -3.0f*(flxY_mid[left]/m_left) +4.0f*(flxY_mid[east]/m_east)-1.0f*(flxY_mid[east]/m_easteast) )/(2.0f*dx); //forward difference
		velY_dy_mid[left] = ( (flxY_mid[north]/m_north)-(flxY_mid[south]/m_south) )/(2.0f*dy); //OK

	}
	for(int j=1; j<=(Ny-1)-2;j++){ //lado direito da rede principal ou seja i=((Nx-1)-1)
		float m_rigth, m_westwest;
		int right = ((Nx-1)-1) + j*(Nx-1);
		int westwest = right-2;
		west=right-1;
		north=right+(Nx-1);
		south=right-(Nx-1);
		m_north = DensityToMass(den_mid[north]);
		m_south = DensityToMass(den_mid[south]);
		m_rigth= DensityToMass(den_mid[right]);
		m_west = DensityToMass(den_mid[west]);
		m_westwest = DensityToMass(den_mid[westwest]);
		velX_dx_mid[right] = ( 3.0f*(flxX_mid[right]/m_rigth) -4.0f*(flxX_mid[west]/m_west) +1.0f*(flxX_mid[westwest]/m_westwest) )/(2.0f*dx); //backwar difference
		velX_dy_mid[right] = ( (flxX_mid[north]/m_north)-(flxX_mid[south]/m_south) )/(2.0f*dy); //OK
		velY_dx_mid[right] = ( 3.0f*(flxY_mid[right]/m_rigth) -4.0f*(flxY_mid[west]/m_west) +1.0f*(flxY_mid[westwest]/m_westwest) )/(2.0f*dx);
		velY_dy_mid[right] = ( (flxY_mid[north]/m_north)-(flxY_mid[south]/m_south) )/(2.0f*dy); //OK
	}
//os 4 cantos em que ambas a derivadas nao podem ser centradas
	int ks;
// i=0 j=0 forward x forward y
	ks = 0 + 0*(Nx-1);
	float m_0 = DensityToMass(den_mid[ks]);
	float m_x1 = DensityToMass(den_mid[ks+1]);
	float m_x2 = DensityToMass(den_mid[ks+2]);
	float m_y1 = DensityToMass(den_mid[ks+(Nx-1)]);
	float m_y2 = DensityToMass(den_mid[ks+2*(Nx-1)]);
	velX_dx_mid[ks] = ( -3.0f*(flxX_mid[ks]/m_0 ) +4.0f*(flxX_mid[ks+1]/m_x1 )-1.0f*(flxX_mid[ks+2]/m_x2 )  )/(2.0f*dx);
	velX_dy_mid[ks] = ( -3.0f*(flxX_mid[ks]/m_0 ) +4.0f*(flxX_mid[ks+1*(Nx-1)]/m_y1 )-1.0f*(flxX_mid[ks+2*(Nx-1)]/m_y2 )  )/(2.0f*dy);
	velY_dx_mid[ks] = ( -3.0f*(flxY_mid[ks]/m_0 ) +4.0f*(flxY_mid[ks+1]/m_x1 )-1.0f*(flxY_mid[ks+2]/m_x2 )  )/(2.0f*dx);
	velY_dy_mid[ks] = ( -3.0f*(flxY_mid[ks]/m_0 ) +4.0f*(flxY_mid[ks+1*(Nx-1)]/m_y1 )-1.0f*(flxY_mid[ks+2*(Nx-1)]/m_y2 )  )/(2.0f*dy);

// i=((Nx-1)-1) j=0 backward x forward y
	ks = ((Nx-1)-1) + 0*(Nx-1);
	m_0 = DensityToMass(den_mid[ks]);
	m_x1 = DensityToMass(den_mid[ks-1]);
	m_x2 = DensityToMass(den_mid[ks-2]);
	m_y1 = DensityToMass(den_mid[ks+(Nx-1)]);
	m_y2 = DensityToMass(den_mid[ks+2*(Nx-1)]);
	velX_dx_mid[ks] = (  3.0f*(flxX_mid[ks]/m_0 ) -4.0f*(flxX_mid[ks-1]/m_x1 )+1.0f*(flxX_mid[ks-2]/m_x2 )  )/(2.0f*dx);
	velX_dy_mid[ks] = ( -3.0f*(flxX_mid[ks]/m_0 ) +4.0f*(flxX_mid[ks+(Nx-1)]/m_y1 )-1.0f*(flxX_mid[ks+2*(Nx-1)]/m_y2 )  )/(2.0f*dy);
	velY_dx_mid[ks] = (  3.0f*(flxY_mid[ks]/m_0 ) -4.0f*(flxY_mid[ks-1]/m_x1 )+1.0f*(flxY_mid[ks-2]/m_x2 )  )/(2.0f*dx);
	velY_dy_mid[ks] = ( -3.0f*(flxY_mid[ks]/m_0 ) +4.0f*(flxY_mid[ks+(Nx-1)]/m_y1 )-1.0f*(flxY_mid[ks+2*(Nx-1)]/m_y2 )  )/(2.0f*dy);

// i=0 j=((Ny-1)-1) forward x backward y
	ks = 0 + ((Ny-1)-1)*(Nx-1);
	m_0 = DensityToMass(den_mid[ks]);
	m_x1 = DensityToMass(den_mid[ks+1]);
	m_x2 = DensityToMass(den_mid[ks+2]);
	m_y1 = DensityToMass(den_mid[ks-(Nx-1)]);
	m_y2 = DensityToMass(den_mid[ks-2*(Nx-1)]);
	velX_dx_mid[ks] = ( -3.0f*(flxX_mid[ks]/m_0 ) +4.0f*(flxX_mid[ks+1]/m_x1 )-1.0f*(flxX_mid[ks+2]/m_x1 )  )/(2.0f*dx);
	velX_dy_mid[ks] = (  3.0f*(flxX_mid[ks]/m_0 ) -4.0f*(flxX_mid[ks-(Nx-1)]/m_y1 )+1.0f*(flxX_mid[ks-2*(Nx-1)]/m_y2 )  )/(2.0f*dy);
	velY_dx_mid[ks] = ( -3.0f*(flxY_mid[ks]/m_0 ) +4.0f*(flxY_mid[ks+1]/m_x1 )-1.0f*(flxY_mid[ks+2]/m_x2 )  )/(2.0f*dx);
	velY_dy_mid[ks] = (  3.0f*(flxY_mid[ks]/m_0 ) -4.0f*(flxY_mid[ks-(Nx-1)]/m_y1 )+1.0f*(flxY_mid[ks-2*(Nx-1)]/m_y2 )  )/(2.0f*dy);

// i=((Nx-1)-1) j=((Ny-1)-1) backward x backward y
	ks = ((Nx-1)-1) + ((Ny-1)-1)*(Nx-1);
	m_0 = DensityToMass(den_mid[ks]);
	m_x1 = DensityToMass(den_mid[ks-1]);
	m_x2 = DensityToMass(den_mid[ks-2]);
	m_y1 = DensityToMass(den_mid[ks-(Nx-1)]);
	m_y2 = DensityToMass(den_mid[ks-2*(Nx-1)]);
	velX_dx_mid[ks] = ( 3.0f*(flxX_mid[ks]/m_0 ) -4.0f*(flxX_mid[ks-1]/m_x1 )+1.0f*(flxX_mid[ks-2]/m_x2 ) )/(2.0f*dx);
	velX_dy_mid[ks] = ( 3.0f*(flxX_mid[ks]/m_0 ) -4.0f*(flxX_mid[ks-(Nx-1)]/m_y1 )+1.0f*(flxX_mid[ks-2*(Nx-1)]/m_y2 ) )/(2.0f*dy);
	velY_dx_mid[ks] = ( 3.0f*(flxY_mid[ks]/m_0 ) -4.0f*(flxY_mid[ks-1]/m_x1 )+1.0f*(flxY_mid[ks-2]/m_x2 ) )/(2.0f*dx);
	velY_dy_mid[ks] = ( 3.0f*(flxY_mid[ks]/m_0 ) -4.0f*(flxY_mid[ks-(Nx-1)]/m_y1 )+1.0f*(flxY_mid[ks-2*(Nx-1)]/m_y2 ) )/(2.0f*dy);

}

float Fluid2D::DensityToMass(float density) {
	return density;
}
