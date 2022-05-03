/************************************************************************************************\
* 2020 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/

#include "includes/Fluid1DLib.h"
#include "includes/Cell1DLib.h"
#include "includes/SetUpParametersLib.h"



using namespace H5;
using namespace std;


Fluid1D::Fluid1D(const SetUpParameters &input_parameters) : TethysBase{input_parameters.SizeX, 0, 1}{
	Nx = input_parameters.SizeX;
	vel_snd = input_parameters.SoundVelocity;
	kin_vis = input_parameters.ShearViscosity;
	col_freq=input_parameters.CollisionFrequency;

	param = {vel_snd,0.0f,0.0f,kin_vis,0.0f,0.0f,col_freq,0.0f};

	char buffer [50];
	sprintf (buffer, "S=%.2fvis=%.2f", vel_snd, kin_vis);
	file_infix = buffer;
	Den = new float[Nx]();
	Vel = new float[Nx]();
	GradVel= new float[Nx]();
	Cur = new float[Nx]();
	vel_snd_arr = new float[Nx]();

	Umain = new StateVec[Nx]();
	Uaux = new StateVec[Nx]();
	Umid = new StateVec[Nx-1]();

}	

Fluid1D::~Fluid1D() = default;


float Fluid1D::VelocityFlux( StateVec U) {
	return 0.5f*U.v()*U.v()+U.n()*U.S()*U.S();
}



float Fluid1D::DensityFlux(StateVec U) {
	return U.n()*U.v();
}

float  Fluid1D::DensitySource( __attribute__((unused)) float n, __attribute__((unused)) float v, __attribute__((unused)) float s){
	return 0;
}
float Fluid1D::VelocitySource(float n, float v, float s, float d3den) {
	return 0;
}

float  Fluid1D::DensitySource(StateVec U){
	return 0;
}
float Fluid1D::VelocitySource(StateVec U) {
	return 0;
}


void Fluid1D::CflCondition(){
		dx = lengX / ( float ) ( Nx - 1 );
		dt = dx/10.0f;
}

void Fluid1D::SetSimulationTime(){
	Tmax=5.0f+0.02f*vel_snd+20.0f/vel_snd;
}
		
void Fluid1D::SetSound(){
	for(int i = 0; i<Nx  ;i++){
		vel_snd_arr[i]= vel_snd;
		Umain[i].S()=vel_snd;
		Uaux[i].S()=vel_snd;
	}
	for(int i = 0; i<Nx-1  ;i++){
		Umid[i].S()=vel_snd;
	}
}
void Fluid1D::SetSound(const std::function<float(float)>& func) {
	for(int i = 0; i<Nx  ;i++){
		vel_snd_arr[i]= func(i*dx);
		Umain[i].S()=func(i*dx);
		Uaux[i].S()=func(i*dx);
	}
	for(int i = 0; i<Nx-1  ;i++){
		Umid[i].S()=func((i+0.5f)*dx);
	}
}



void Fluid1D::InitialCondRand(){
	random_device rd;
	float maxrand;
	maxrand = (float) random_device::max();

	for (int i = 0; i < Nx; i++ ){
		float noise = (float) rd()/ maxrand ;
		Umain[i].n()= 1.0f + 0.0001f * (noise - 0.5f);
		Umain[i].v()= 0.0f;
	}
	this->SetSound();
}

void Fluid1D::InitialCondTest(){
 	for (int i = 0; i < Nx; i++ ){
		Umain[i].n()=1.0;
	    Umain[i].v()=(i>Nx/3 && i<2*Nx/3 ) ? 1.0f : 0.1f;
	}
	this->SetSound();
}
void Fluid1D::InitialCondGeneral(function<float(float)> fden, function<float(float)> fvx) {
	float x;
	for (int i = 0; i < Nx; ++i) {
		x=i*dx;
		Umain[i].n()=fden(x);
		Umain[i].v()=fvx(x);
	}
	this->SetSound();
}



void Fluid1D::CreateFluidFile(){
	std::string previewfile = "preview_1D_" + file_infix + ".dat" ;
	data_preview.open (previewfile);
	data_preview << scientific; 
}

void Fluid1D::WriteFluidFile(float t){
	int pos_end = Nx - 1 ;
	int pos_ini = 0;
	if (!isfinite(Umain[pos_end].n()) || !isfinite(Umain[pos_ini].n()) || !isfinite(Umain[pos_end].v()) ||
	    !isfinite(Umain[pos_ini].v())) {
		cerr << "ERROR: numerical method failed to converge" <<"\nExiting"<< endl;
		exit(EXIT_FAILURE);
	}
//data_preview << t << "\t" << Umain[pos_end/2] << "\n";
	data_preview << t << "\t" << Umain[pos_ini] << "\t" << Umain[pos_end] << "\n";
}


void Fluid1D::Richtmyer() {

	//CalcDensityLaplacian(Umain);
	RichtmyerStep1();
	
	//CalcDensityLaplacian(Umid);
	RichtmyerStep2();
}
void Fluid1D::RichtmyerStep1() {
	for ( int i = 0; i <= Nx - 2; i++ ){
		float den_avg   = 0.5f * (Umain[i+1] + Umain[i] ).n();
		float vel_avg   = 0.5f * (Umain[i+1] + Umain[i] ).v();
		Umid[i].n() = den_avg - 0.5f*(dt/dx)*(DensityFlux(Umain[i+1]) - DensityFlux(Umain[i]))
							+ (0.5f*dt) * DensitySource(0.5f * (Umain[i+1] + Umain[i] )) ;
		Umid[i].v() = vel_avg - 0.5f*(dt/dx)*(VelocityFlux(Umain[i+1]) - VelocityFlux(Umain[i]))
							+ (0.5f*dt) * VelocitySource(0.5f * (Umain[i+1] + Umain[i] )) ;
	}
}
void Fluid1D::RichtmyerStep2() {
	for ( int i = 1; i <= Nx - 2; i++ ){
		StateVec Uold(Umain[i]);
		float den_old = Uold.n();
		float vel_old = Uold.v();
//		float den_old = Umain[i].n();
//		float vel_old = Umain[i].v();
		Umain[i].n() = den_old - (dt/dx)*(DensityFlux(Umid[i]) - DensityFlux(Umid[i-1]))
				                 + dt * DensitySource(Uold);
		Umain[i].v() = vel_old - (dt/dx)*(VelocityFlux(Umid[i]) - VelocityFlux(Umid[i-1]))
		                         + dt * VelocitySource(Uold);
	}
}




void Fluid1D::VelocityToCurrent() {
	for(int i=0; i <= Nx - 1; i++){
		Cur[i] = Vel[i] * Den[i];
	}
}

int Fluid1D::GetSnapshotStep() const { return snapshot_step;}
int Fluid1D::GetSnapshotFreq() const {return snapshot_per_period;}


bool Fluid1D::Snapshot() const {
	return TimeStepCounter % snapshot_step == 0;
}


void Fluid1D::SaveSnapShot(){

	hsize_t dim_atr[1] = { 1 };
	DataSpace atr_dataspace = DataSpace (1, dim_atr );

	int points_per_period = static_cast<int>((2.0 * MAT_PI / RealFreq()) / dt);
	snapshot_step = 1; //points_per_period / snapshot_per_period;

	string str_time = to_string(TimeStepCounter );/// snapshot_step);
	str_time.insert(str_time.begin(), 5 - str_time.length(), '0');
	string name_dataset = "snapshot_" + str_time;

	float currenttime= static_cast<float>(TimeStepCounter) * dt;
	CopyFields();

	DataSet dataset_den = GrpDen->createDataSet(name_dataset, HDF5FLOAT, *DataspaceDen);
	Attribute atr_step_den = dataset_den.createAttribute("time step", HDF5INT, atr_dataspace);
	Attribute atr_time_den = dataset_den.createAttribute("time", HDF5FLOAT, atr_dataspace);
	dataset_den.write(Den, HDF5FLOAT);
	dataset_den.close();
	atr_step_den.write(HDF5INT, &TimeStepCounter);
	atr_time_den.write(HDF5FLOAT , &currenttime);
	atr_step_den.close();
	atr_time_den.close();

	DataSet dataset_vel_x = GrpVelX->createDataSet(name_dataset, HDF5FLOAT, *DataspaceVelX);
	Attribute atr_step_vel_x = dataset_vel_x.createAttribute("time step", HDF5INT, atr_dataspace);
	Attribute atr_time_vel_x = dataset_vel_x.createAttribute("time", HDF5FLOAT, atr_dataspace);
	dataset_vel_x.write(Vel, HDF5FLOAT);
	dataset_vel_x.close();
	atr_step_vel_x.write(HDF5INT, &TimeStepCounter);
	atr_time_vel_x.write(HDF5FLOAT , &currenttime);
	atr_step_vel_x.close();
	atr_time_vel_x.close();
}

void Fluid1D::RungeKuttaTVD() {

	float DenNumFluxW;
	float DenNumFluxE;
	float VelNumFluxW;
	float VelNumFluxE;

	float DenNumSourceW;
	float DenNumSourceE;
	float VelNumSourceW;
	float VelNumSourceE;

	StateVec UEleft{};
	StateVec UEright{};
	StateVec UWleft{};
	StateVec UWright{};

	StateVec U_1[Nx];
	StateVec U_2[Nx];

	//CalcDensityLaplacian(Umain,Nx);
	for (int i = 1; i < Nx-1; ++i) {

		// reconstruction process
		CellHandler1D cell(i, this, Umain);
		UEleft  = cell.WENO3(Nx,'E','L');
		UEright = cell.WENO3(Nx,'E','R');
		UWleft  = cell.WENO3(Nx,'W','L');
		UWright = cell.WENO3(Nx,'W','R');

		// calculates flux terms
		DenNumFluxE = NumericalFlux::Central(this,UEleft,UEright).n();
		DenNumFluxW = NumericalFlux::Central(this,UWleft,UWright).n();
		VelNumFluxE = NumericalFlux::Central(this,UEleft,UEright).v();
		VelNumFluxW = NumericalFlux::Central(this,UWleft,UWright).v();

		// calculates source terms
		DenNumSourceE = NumericalSource::Average(this,UEleft,UEright).n();
		DenNumSourceW = NumericalSource::Average(this,UWleft,UWright).n();
		VelNumSourceE = NumericalSource::Average(this,UEleft,UEright).v();
		VelNumSourceW = NumericalSource::Average(this,UWleft,UWright).v();

		float Ln = -(DenNumFluxE-DenNumFluxW)/dx + 0.5f*(DenNumSourceE+DenNumSourceW);
		float Lv = -(VelNumFluxE-VelNumFluxW)/dx + 0.5f*(VelNumSourceE+VelNumSourceW);

		// RK
		U_1[i].n() = Umain[i].n() + dt*Ln;
		U_1[i].v() = Umain[i].v() + dt*Lv;
	}
	// fixes boundary cells for the first iteration (this solution only works with periodic boundary condition)
	U_1[0] = U_1[Nx-2];
	U_1[Nx-1] = U_1[1];

	//CalcDensityLaplacian(U_1,Nx);
	for (int i = 1; i < Nx-1; ++i) {

		// reconstruction process
		CellHandler1D cell(i, this, U_1);
		UEleft  = cell.WENO3(Nx,'E','L');
		UEright = cell.WENO3(Nx,'E','R');
		UWleft  = cell.WENO3(Nx,'W','L');
		UWright = cell.WENO3(Nx,'W','R');

		// calculates flux terms
		DenNumFluxE = NumericalFlux::Central(this,UEleft,UEright).n();
		DenNumFluxW = NumericalFlux::Central(this,UWleft,UWright).n();
		VelNumFluxE = NumericalFlux::Central(this,UEleft,UEright).v();
		VelNumFluxW = NumericalFlux::Central(this,UWleft,UWright).v();

		// calculates source terms
		DenNumSourceE = NumericalSource::Average(this,UEleft,UEright).n();
		DenNumSourceW = NumericalSource::Average(this,UWleft,UWright).n();
		VelNumSourceE = NumericalSource::Average(this,UEleft,UEright).v();
		VelNumSourceW = NumericalSource::Average(this,UWleft,UWright).v();

		float Ln = -(DenNumFluxE-DenNumFluxW)/dx + 0.5f*(DenNumSourceE+DenNumSourceW);
		float Lv = -(VelNumFluxE-VelNumFluxW)/dx + 0.5f*(VelNumSourceE+VelNumSourceW);

		// RK
		U_2[i].n() = 0.75f*Umain[i].n() + 0.25f*U_1[i].n() + 0.25f*dt*Ln;
		U_2[i].v() = 0.75f*Umain[i].v() + 0.25f*U_1[i].v() + 0.25f*dt*Lv;
	}
	// fixes boundary cells for the first iteration (this solution only works with periodic boundary condition)
	U_2[0] = U_2[Nx-2];
	U_2[Nx-1] = U_2[1];

	//CalcDensityLaplacian(U_2,Nx);
	for (int i = 1; i < Nx-1; ++i) {

		// reconstruction process
		CellHandler1D cell(i, this, U_2);
		UEleft  = cell.WENO3(Nx,'E','L');
		UEright = cell.WENO3(Nx,'E','R');
		UWleft  = cell.WENO3(Nx,'W','L');
		UWright = cell.WENO3(Nx,'W','R');

		// calculates flux terms
		DenNumFluxE = NumericalFlux::Central(this,UEleft,UEright).n();
		DenNumFluxW = NumericalFlux::Central(this,UWleft,UWright).n();
		VelNumFluxE = NumericalFlux::Central(this,UEleft,UEright).v();
		VelNumFluxW = NumericalFlux::Central(this,UWleft,UWright).v();

		// calculates source terms
		DenNumSourceE = NumericalSource::Average(this,UEleft,UEright).n();
		DenNumSourceW = NumericalSource::Average(this,UWleft,UWright).n();
		VelNumSourceE = NumericalSource::Average(this,UEleft,UEright).v();
		VelNumSourceW = NumericalSource::Average(this,UWleft,UWright).v();

		float Ln = -(DenNumFluxE-DenNumFluxW)/dx + 0.5f*(DenNumSourceE+DenNumSourceW);
		float Lv = -(VelNumFluxE-VelNumFluxW)/dx + 0.5f*(VelNumSourceE+VelNumSourceW);

		// RK
		Umain[i].n() = (1.0f/3.0f)*Umain[i].n() + (2.0f/3.0f)*U_2[i].n() + (2.0f/3.0f)*dt*Ln;
		Umain[i].v() = (1.0f/3.0f)*Umain[i].v() + (2.0f/3.0f)*U_2[i].v() + (2.0f/3.0f)*dt*Lv;
	}
}

void Fluid1D::McCormack() {
	for (int i = 1; i < Nx-1; ++i) {
		Uaux[i].n()=Umain[i].n()-(dt/dx)*(DensityFlux(Umain[i+1])-DensityFlux(Umain[i]));
		Uaux[i].v()=Umain[i].v()-(dt/dx)*(VelocityFlux(Umain[i+1])-VelocityFlux(Umain[i]));
	}
	for (int i = 1; i < Nx-1; ++i) {
		Umain[i].n()=0.5f*(Umain[i].n()+Uaux[i].n())-(0.5f*dt/dx)*(DensityFlux(Uaux[i])-DensityFlux(Uaux[i-1]));
		Umain[i].v()=0.5f*(Umain[i].v()+Uaux[i].v())-(0.5f*dt/dx)*(VelocityFlux(Uaux[i])-VelocityFlux(Uaux[i-1]));
	}
}

void Fluid1D::Upwind(){
	for (int i = 1; i < Nx-1; ++i) {
		Umain[i].n()=Umain[i].n()-(dt/dx)*(DensityFlux(Umain[i])-DensityFlux(Umain[i-1]));
		Umain[i].v()=Umain[i].v()-(dt/dx)*(VelocityFlux(Umain[i])-VelocityFlux(Umain[i-1]));
	}
}

void Fluid1D::LaxFriedrichs(){
	for (int i = 1; i < Nx-1; ++i) {
		Umain[i].n()=Den[i];
		Umain[i].v()=Vel[i];
		Den[i]=0.5f*(Umain[i-1].n()+Umain[i+1].n())-0.5f*(dt/dx)*(DensityFlux(Umain[i+1])-DensityFlux(Umain[i-1]));
		Vel[i]=0.5f*(Umain[i-1].v()+Umain[i+1].v())-0.5f*(dt/dx)*(VelocityFlux(Umain[i+1])-VelocityFlux(Umain[i-1]));
	}
}

float Fluid1D::JacobianSpectralRadius(StateVec U) {
	float l1=abs(U.v()+ vel_snd*sqrt(U.n()));
	float l2=abs(U.v()- vel_snd*sqrt(U.n()));
	return max(l1,l2);
}

StateVec Fluid1D::ConservedFlux(StateVec U) {
	StateVec Uout{};
	Uout.n()= this->DensityFlux(U);
	Uout.v()= this->VelocityFlux(U);
	return Uout;
}

float Fluid1D::JacobianSignum( StateVec U, std::string key) {

	float l1= Signum(U.v()+vel_snd*sqrt(U.n()));
	float l2= Signum(U.v()-vel_snd*sqrt(U.n()));
	float entry=0;
	if(key=="11"){
		entry=l1*0.5f;
	}else if(key=="12"){
		entry=l2* sqrt(U.n())/(2.0f*vel_snd);
	}else if(key=="21"){
		entry=l1* vel_snd/sqrt(U.n());
	}else if(key=="22"){
		entry=l2*0.5f;
	}else entry=0.0f;
return entry;
}

void Fluid1D::CopyFields() {
	for (int i = 0; i < Nx; ++i) {
		Den[i]=Umain[i].n();
		Vel[i]=Umain[i].v();
		Cur[i]=Umain[i].v()*Umain[i].n();
	}
}


void Fluid1D::SaveSound() {
	DataSet dataset_vel_snd = GrpDat->createDataSet("Sound velocity", HDF5FLOAT, *DataspaceVelSnd);
	dataset_vel_snd.write(vel_snd_arr, HDF5FLOAT);
	dataset_vel_snd.close();
}


void Fluid1D::CalcDensityLaplacian(StateVec* u_vec) {

	// calculates the laplacian at the extreme cells
	u_vec[0].lap_n()    = (2.0f*u_vec[0].n() - 5.0f*u_vec[1].n() + 4.0f*u_vec[2].n() - u_vec[2].n()) / (dx*dx);
	u_vec[Nx-1].lap_n() = (2.0f*u_vec[Nx-1].n() - 5.0f*u_vec[Nx-2].n() + 4.0f*u_vec[Nx-3].n() - u_vec[Nx-4].n()) / (dx*dx);

	// calculates the laplacian for all the other cells
	for(int i = 1; i < Nx-1; ++i) {
		u_vec[i].lap_n() = (u_vec[i-1].n() - 2.0f*u_vec[i].n() + u_vec[i+1].n()) / (dx*dx);
	}
}