#include "TethysLib.h"


using namespace H5;
using namespace std;





/*....................................................................*/
/*........ General Functions .........................................*/
/*....................................................................*/
void TethysBase::BannerDisplay() const {
cout<<"\n" ;
	cout<<"╔═════════════════════════════════════════════════════════════════════════╗\n";
	cout<<"║\033[2m  ▆▆▆▆▆▆▆▆▆▆▆ ▆▆▆▆▆▆▆▆▆▆  ▆▆▆▆▆▆▆▆▆▆▆ ▆▆▆▆▆ ▆▆▆▆▆ ▆▆▆▖   ▗▆▆▆ ▗▆▆▆▆▆▆▆▖  \033[0m║\n";
	cout<<"║\033[2m  █▘  ▐█▌  ▝█  ▐█▌    ▝█  █▘  ▐█▌  ▝█  ▐█▌   ▐█▌     █▌ ▐█   ▐█▌     ▝█  \033[0m║\n";
	cout<<"║\033[2m      ▐█▌      ▐█▌▆▆▆█        ▐█▌      ▐█▌▆▆▆▐█▌      ▐█▌     ▝██▆▆▆▆▆▖  \033[0m║\n";
	cout<<"║\033[2m      ▐█▌      ▐█▌    ▗▉      ▐█▌      ▐█▌   ▐█▌      ▐█▌    ▗       ██  \033[0m║\n";
	cout<<"║\033[2m     ▆███▆    ▆███▆▆▆██▉     ▆███▆    ▆███▆ ▆███▆    ▆███▆   ▐█▆▆▆▆▆██▘  \033[0m║\n";
	cout<<"║                                                                         ║\n";
	cout<<"║ \033[1mTwo-dimensional Emitter of THz, Hydrodynamic Simulation.  Version 2.2.1\033[0m ║\n";
	cout<<"╚═════════════════════════════════════════════════════════════════════════╝\n";
}

void TethysBase::WelcomeScreen() const {
	cout << "\nFermi velocity\t\033[1mvF\t" << vel_fer << " v\342\202\200\033[0m\n";
	if (this->PhaseVel() < vel_fer) {
		cout << "Phase velocity\t\033[1mS'\t" << this->PhaseVel()
		     << " v\342\202\200\033[0m  \033[1;5;7;31m WARNING plasmon in damping region \033[0m" << endl;
	} else {
		cout << "Phase velocity\t\033[1mS'\t" << this->PhaseVel() << " v\342\202\200\033[0m\n";
	}
	cout << "Viscosity \t\033[1m\316\267\t" << kin_vis << "\033[0m\n";
	if (kin_vis != 0.0) {
		if (2.0f * kin_vis * dt > 0.95 * dx * dx * dy * dy / (dx * dx + dy * dy)) {
			cout << "Reynolds n. \t\033[1mRe\t" << 1.0 / kin_vis
			     << "\033[0m \033[1;5;7;31m WARNING ftcs scheme may not converge.  \033[0m" << endl;
		} else {
			cout << "Reynolds n. \t\033[1mRe\t" << 1.0 / kin_vis << "\033[0m\n";
		}
	}
	cout << "Collision rate \t\033[1m\316\275\t" << col_freq << " v\342\202\200/L\033[0m\n";
	cout << "Cyclotron frequency \t\033[1m\317\211c\t" << cyc_freq << " v\342\202\200/L\n\033[0m\n";
	cout << "Theoretical frequency \033[1m\317\211=\317\211'+i\317\211''\033[0m\n";
	cout << "\033[1m\317\211'\t" << this->RealFreq() << " v\342\202\200/L\t2\317\200/\317\211'\t" << 2.0 * MAT_PI /
	                                                                                                 this->RealFreq()
	     << " L/v\342\202\200\033[0m\n";
	cout << "\033[1m\317\211''\t" << this->ImagFreq() << " v\342\202\200/L\t2\317\200/\317\211''\t" << 2.0 * MAT_PI /
	                                                                                                   this->ImagFreq()
	     << " L/v\342\202\200\033[0m\n";
	cout << "\nDetermined maximum simulated time\t\033[1m\nT\342\202\230\342\202\220\342\202\223\t" << Tmax
	     << " L/v\342\202\200\t\342\211\210" << Tmax / dt << "\033[0m\t time steps" << endl;
	cout << "Discretisation\n";
	if (2==RANK){
	cout << "\033[1m\316\224t\t" << dt << " L/v\342\202\200\t\316\224x\t" << dx << " L\t" << "\316\224y\t" << dy
	     << " L\033[0m\n";
	cout << "Simulation grid\t\033[1m" << Nx - 1 << " x " << Ny - 1 << "\033[0m\n";
	}
	if (1==RANK){
		cout << "\033[1m\316\224t\t" << dt << " L/v\342\202\200\t\316\224x\t" << dx << " L\033[0m\n";
		cout << "Simulation grid\t\033[1m" << Nx - 1 << "\033[0m\n";
	}
}



std::string TethysBase::GetInfix() const {return file_infix;}
float TethysBase::GetTmax() const{return Tmax;}
int TethysBase::Rank() const{ return RANK; }
int TethysBase::SizeX() const{ return Nx; }
int TethysBase::SizeY() const{ return Ny; }
float TethysBase::GetVelSnd() const{ return vel_snd; }
float TethysBase::GetKinVis() const{ return kin_vis; }
float TethysBase::GetColFreq() const{ return col_freq; }
float TethysBase::GetVelFer() const{ return vel_fer;  }
float TethysBase::GetCycFreq() const{ return cyc_freq; }
float TethysBase::GetDx() const{return dx;}
float TethysBase::GetDy() const{return dy;}
float TethysBase::GetDt() const{return dt;}
float TethysBase::GetLengthX() const{return lengX;}
float TethysBase::GetLengthY() const{return lengY;}



void TethysBase::SetTmax(float x){ Tmax=x;}
void TethysBase::SetVelSnd(float x){ vel_snd=x; }
void TethysBase::SetKinVis(float x){ kin_vis=x;}
void TethysBase::SetColFreq(float x){ col_freq=x; }
void TethysBase::SetVelFer(float x){ vel_fer=x;}
void TethysBase::SetCycFreq(float x) { cyc_freq=x;}
void TethysBase::SetDx(float x){ dx=x;}
void TethysBase::SetDy(float x){ dy=x;}
void TethysBase::SetDt(float x){ dt=x;}
void TethysBase::SetLengthX(float x){lengX=x;}
void TethysBase::SetLengthY(float x){lengY=x;}


void TethysBase::CloseHdf5File() const{
	DataspaceVelX->close();
	DataspaceDen->close();
	DataspaceVelSnd->close();
	GrpDat->close();
	GrpDen->close();
	GrpVelX->close();
	if (RANK == 2) {
		GrpVelY->close();
		DataspaceVelY->close();
	}
	GrpDat->close();
	Hdf5File->close();
}


void TethysBase::WriteAttributes(){
	int total_steps= static_cast<int>(Tmax / dt);
	//Create the data space for the attribute.
	hsize_t dim_atr[1] = { 1 };
	DataSpace atr_dataspace = DataSpace (1, dim_atr );
	// Create a group attribute.
	Attribute atr_vel_snd  = GrpDat->createAttribute("Sound velocity", HDF5FLOAT, atr_dataspace);
	Attribute atr_kin_vis = GrpDat->createAttribute("Kinetic viscosity", HDF5FLOAT, atr_dataspace);
	Attribute atr_col_freq = GrpDat->createAttribute("Collision frequency", HDF5FLOAT, atr_dataspace);
	Attribute atr_cyc_freq  = GrpDat->createAttribute("Cyclotron frequency", HDF5FLOAT, atr_dataspace);
	Attribute atr_vel_fer  = GrpDat->createAttribute("Fermi velocity", HDF5FLOAT, atr_dataspace);
	Attribute atr_dx = GrpDat->createAttribute("Space discretisation step x", HDF5FLOAT, atr_dataspace);
	Attribute atr_dt = GrpDat->createAttribute("Time discretisation step", HDF5FLOAT, atr_dataspace);
	Attribute atr_total_time = GrpDat->createAttribute("Total simulation time", HDF5FLOAT, atr_dataspace);
	Attribute atr_num_space_points_x = GrpDat->createAttribute("Number of spatial points x", HDF5INT, atr_dataspace);
	Attribute atr_num_time_steps = GrpDat->createAttribute("Number of time steps", HDF5INT, atr_dataspace);
	Attribute atr_dy;
	Attribute atr_num_space_points_y;
	Attribute atr_aspect_ratio;
	if (RANK == 2) {
		atr_aspect_ratio = GrpDat->createAttribute("Aspect ratio", HDF5FLOAT, atr_dataspace);
		atr_dy = GrpDat->createAttribute("Space discretisation step y", HDF5FLOAT, atr_dataspace);
		atr_num_space_points_y = GrpDat->createAttribute("Number of spatial points y", HDF5INT,
		                                                 atr_dataspace);
	}
	// Write the attribute data.
	atr_vel_snd.write(HDF5FLOAT, &vel_snd);
	atr_vel_fer.write(HDF5FLOAT, &vel_fer);
	atr_cyc_freq.write(HDF5FLOAT, &cyc_freq);
	atr_col_freq.write(HDF5FLOAT, &col_freq);
	atr_kin_vis.write(HDF5FLOAT, &kin_vis);
	atr_dx.write(HDF5FLOAT, &dx);
	atr_dt.write(HDF5FLOAT, &dt);
	atr_num_space_points_x.write(HDF5INT, &Nx);
	atr_total_time.write(HDF5FLOAT, &Tmax);
	atr_num_time_steps.write(HDF5INT, &total_steps);
	if (RANK == 2) {
		float aspect_ratio= lengX/lengY;
		atr_dy.write(HDF5FLOAT, &dy);
		atr_num_space_points_y.write(HDF5INT, &Ny);
		atr_aspect_ratio.write(HDF5FLOAT,&aspect_ratio);
	}
	// Close the attributes.
	atr_num_time_steps.close();
	atr_col_freq.close();
	//atr_vel_fer.close();
	atr_vel_snd.close();
	atr_kin_vis.close();
	atr_dx.close();
	atr_dt.close();
	atr_total_time.close();
	atr_num_space_points_x.close();
	if (RANK == 2) {
		atr_num_space_points_y.close();
		atr_dy.close();
	}
}


TethysBase::~TethysBase(){
	if(HDF5fileOpen) {
		this->CloseHdf5File();
		delete GrpDat;
		delete GrpDen;
		delete GrpVelX;
		delete DataspaceDen;
		delete DataspaceVelX;
		delete DataspaceVelSndMid;
		delete DataspaceVelSnd;
		if (RANK == 2) {
			delete GrpVelY;
			delete DataspaceVelY;
		}
		delete Hdf5File;
	}
}







bool TethysBase::HDF5fileOpen=false;
int TethysBase::TimeStepCounter=0;
float TethysBase::TimeStamp=0.0f;

TethysBase::TethysBase(int size_nx, int size_ny, int dimension){
	Nx = size_nx;
	Ny = size_ny;
	RANK=dimension;

	char buffer [50];
	if(RANK==1) {
		sprintf(buffer, "Fluido1D_Nx=%d", Nx);
	}
	if(RANK==2) {
		sprintf(buffer, "Fluido2D_Nx=%d_Ny=%d", Nx, Ny);
	}
	file_infix = buffer;

	if(RANK==1){
		hsize_t dimsf[1];// dataset dimensions
		dimsf[0] = static_cast<hsize_t>(Nx);
		DataspaceDen = new DataSpace(RANK, dimsf );
		DataspaceVelX = new DataSpace(RANK, dimsf );
		DataspaceVelSnd = new DataSpace(RANK, dimsf );
	}
	if(RANK==2){
		hsize_t dimsf[2];
		dimsf[0] = static_cast<hsize_t>(Ny);
		dimsf[1] = static_cast<hsize_t>(Nx);  //troquei !
		DataspaceVelSnd = new DataSpace(RANK, dimsf );
		DataspaceDen = new DataSpace(RANK, dimsf );
		DataspaceVelX = new DataSpace(RANK, dimsf );
		DataspaceVelY = new DataSpace(RANK, dimsf );
		dimsf[0] = static_cast<hsize_t>(Ny-1);
		dimsf[1] = static_cast<hsize_t>(Nx-1);  //troquei !
		DataspaceVelSndMid = new DataSpace(RANK, dimsf );
	}

	Hdf5File = nullptr;
	GrpDat = nullptr;
	GrpDen = nullptr;
	GrpVelX = nullptr;
	GrpVelY = nullptr;
}


void TethysBase::CreateHdf5File(){
	std::string hdf5name;
	TethysBase::HDF5fileOpen=true;
	if(RANK==1){
		hdf5name = "hdf5_1D_" + this->GetInfix() + ".h5" ;
	}
	if(RANK==2){
		hdf5name = "hdf5_2D_" + this->GetInfix() + ".h5" ;
	}
	H5std_string  file_name(hdf5name );
	Hdf5File = new H5File(file_name, H5F_ACC_TRUNC );
	GrpDat = new Group(Hdf5File->createGroup("/Data" ));
	GrpDen = new Group(Hdf5File->createGroup("/Data/Density" ));
	GrpVelX = new Group(Hdf5File->createGroup("/Data/VelocityX" ));
	if(RANK==2) {
		GrpVelY = new Group(Hdf5File->createGroup("/Data/VelocityY"));
	}
}


void TethysBase::OpenHdf5File(const std::string& hdf5name){
	TethysBase::HDF5fileOpen=true;
	Hdf5File = new H5File(hdf5name, H5F_ACC_RDONLY );
	GrpDat = new Group(Hdf5File->openGroup("/Data" ));
	GrpDen = new Group(Hdf5File->openGroup("/Data/Density" ));
	GrpVelX = new Group(Hdf5File->openGroup("/Data/VelocityX" ));
	if(RANK==2) {
		GrpVelY = new Group(Hdf5File->openGroup("/Data/VelocityY"));
	}
}




///////////////////////////////////////////////////////////////////////////////////////////////////////77
/*void Record_Log_File(float vel_snd, float vel_fer, float col_freq, float dt, float dx, float dy, float tmax){
	ofstream logfile;
	logfile.open("Simulation.log",std::ios_base::app);
	time_t time_raw;
	struct tm * time_info;
	time (&time_raw);
	time_info = localtime (&time_raw);
	char buffer [80];
	strftime (buffer,80,"%F %H:%M:%S\n",time_info);
	logfile << "\n#Simulation @ " << buffer ;
	logfile << "#parameters:\n";
	logfile << "#vel_snd \t vel_fer \t col_freq  \t w' \t w'' \n";
	logfile << vel_snd << "\t" << vel_fer << "\t" << col_freq << "\t" << RealFreq(vel_snd, vel_fer, col_freq, 1) << "\t" << ImagFreq(
			vel_snd, vel_fer, col_freq) << "\n";
	logfile << "#discretisation:\n";
	logfile << "#dt\tdx\ttmax\ttime steps\tspace points\n";
	logfile << dt << "\t" << dx << "\t" << dy << "\t" << tmax << "\t" << (int) (tmax / dt) << "\t" << (int) 1 / dx << endl;
}
*/

float TethysBase::PhaseVel()const{
	return sqrt(vel_snd*vel_snd+0.5f*vel_fer*vel_fer + 0.0625f );
}
float TethysBase::RealFreq()const {
	float mode = 1.0f;
	float vel_phs = this->PhaseVel();
	float vel_phs_sqr = vel_phs*vel_phs ;
	if (1.0f < vel_phs ){
		mode = 2.0f*mode-1.0f;
	}
	else{
		mode = 2.0f*mode;
		}
	return fabs(vel_phs_sqr - 0.5625f ) * MAT_PI * mode / (2.0f * vel_phs );
}
float TethysBase::ImagFreq()const {
	float vel_phs = this->PhaseVel();
	float vel_phs_sqr = vel_phs*vel_phs ;
	return (vel_phs_sqr - 0.5625f ) * log(fabs( (vel_phs+0.75f)/(vel_phs-0.75f) )) / (2.0f * vel_phs ) - col_freq*(1.0f-0.125f/vel_phs);
}


SetUpParameters::SetUpParameters() {
	SizeX=101;
	SizeY=101;
	SoundVelocity = 30.0f;
	FermiVelocity = 10.0f;
	CollisionFrequency = 0.01f;
	ShearViscosity = 0.0f;
	CyclotronFrequency = 0.0f;
	SaveMode = 1;
	this->DefineGeometry();
}


SetUpParameters::SetUpParameters(float sound, float fermi, float coll, float visco, float cyclo, int mode, float aspect){
	SizeX=101;
	SizeY=101;
	try {
		SoundVelocity = sound;
		FermiVelocity = fermi;
		CollisionFrequency = coll;
		ShearViscosity = visco;
		CyclotronFrequency = cyclo;
		SaveMode = mode;
		AspectRatio = aspect;
		this->ExceptionsChecking();
	}catch (const char* msg) {
		cerr << msg <<"\nExiting"<< endl;
		exit(EXIT_FAILURE);
	}
	this->DefineGeometry();
}


SetUpParameters::SetUpParameters(int argc, char ** argv) {
	SizeX=201;
	SizeY=201;
	if(argc==7||argc==8){
		try {
			SoundVelocity = strtof(argv[1], nullptr);
			FermiVelocity = strtof(argv[2], nullptr);
			CollisionFrequency = strtof(argv[3], nullptr);
			ShearViscosity = strtof(argv[4], nullptr);
			CyclotronFrequency = strtof(argv[5], nullptr);
			SaveMode = (int) strtol(argv[6], nullptr, 10);    // full data or light save option
			if (argc == 8) {
				AspectRatio = strtof(argv[7], nullptr);
			}
			this->ExceptionsChecking();
		}catch (const char* msg) {
			cerr << msg <<"\nExiting"<< endl;
			exit(EXIT_FAILURE);
		}
	}
	else{
		try {
			cout << "Define S value: "; // throw exceptions if the velocities or frequency are negative or if S<Vf
			cin >> SoundVelocity;
			cout << "Define vF value: ";
			cin >> FermiVelocity;
			cout << "Define kinetic viscosity: ";
			cin >> ShearViscosity;
			cout << "Define collision frequency: ";
			cin >> CollisionFrequency;
			cout << "Define cyclotron frequency: ";
			cin >> CyclotronFrequency;
			cout << "Define the aspect ratio x:y ";
			cin >> AspectRatio;
			cout << "Define data_save_mode value (0-> light save | 1-> full data): ";
			cin >> SaveMode;
			this->ExceptionsChecking();
		}catch (const char* msg) {
			cerr << msg  <<"\nExiting"<< endl;
			exit(EXIT_FAILURE);
		}
	}
	this->DefineGeometry();
}



void SetUpParameters::ExceptionsChecking() const{
	if(SoundVelocity<=0.0f){
		throw "ERROR: Unphysical Sound Velocity";
	}
	if(FermiVelocity<=0.0f){
		throw "ERROR: Unphysical Fermi Velocity";
	}
	if(ShearViscosity<0.0f){
		throw "ERROR: Unphysical Shear Viscosity";
	}
	if(CollisionFrequency<0.0f){
		throw "ERROR: Unphysical Collision Frequency";
	}
	if(CyclotronFrequency<0.0f){
		throw "ERROR: Unphysical Cyclotron Frequency";
	}
	if( SaveMode != 0 && SaveMode != 1  ) {
		throw "ERROR: Unknown save mode option";
	}
}

void SetUpParameters::ParametersFromHdf5File(const std::string& hdf5name){
	H5File *hdf5_file;
	Group *grp_dat;
	try{
		Exception::dontPrint();
		hdf5_file = new H5File(hdf5name, H5F_ACC_RDONLY);
		grp_dat = new Group(hdf5_file->openGroup("/Data"));
	}
	catch( FileIException &file_error )
	{
		cerr<<"Unable to open HDF5 file\t"<< file_error.getDetailMsg() <<"\nExiting\n";
		exit(EXIT_FAILURE);
	}
	catch (...) {
		cerr<<"Unknown error found\nExiting";
		exit(EXIT_FAILURE);
	}
	Attribute attr_n_x(grp_dat->openAttribute("Number of spatial points x"));
	attr_n_x.read(attr_n_x.getDataType(), &SizeX);
	Attribute attr_n_y(grp_dat->openAttribute("Number of spatial points y"));
	attr_n_y.read(attr_n_y.getDataType(), &SizeY);
	Attribute attr_snd(grp_dat->openAttribute("Sound velocity"));
	attr_snd.read(attr_snd.getDataType(), &SoundVelocity);
	Attribute attr_fer(grp_dat->openAttribute("Fermi velocity"));
	attr_fer.read(attr_fer.getDataType(), &FermiVelocity);
	Attribute attr_vis(grp_dat->openAttribute("Kinetic viscosity"));
	attr_vis.read(attr_vis.getDataType(), &ShearViscosity);
	Attribute attr_cyc(grp_dat->openAttribute("Cyclotron frequency"));
	attr_cyc.read(attr_cyc.getDataType(), &CyclotronFrequency);
	Attribute attr_col(grp_dat->openAttribute("Collision frequency"));
	attr_col.read(attr_col.getDataType(), &CollisionFrequency);
	AspectRatio = (static_cast<float>(SizeX-1)) / (static_cast<float>(SizeY-1));
	attr_n_x.close();
	attr_n_y.close();
	attr_snd.close();
	attr_fer.close();
	attr_vis.close();
	attr_cyc.close();
	attr_col.close();
	grp_dat->close();
	hdf5_file->close();
}

void SetUpParameters::DefineGeometry() {
	if(AspectRatio>1.0f){
		Length=1.0f*AspectRatio;
		Width=1.0f;
		SizeY=201;
		SizeX= static_cast<int>( static_cast<float>(SizeY-1)*AspectRatio)+1;
	}
	if(AspectRatio==1.0f){
		Length=1.0f;
		Width=1.0f;
		SizeX=201;
		SizeY=201;
	}
	if(AspectRatio<1.0f){
		Length=1.0f;
		Width=1.0f/AspectRatio;
		SizeX=201;
		SizeY= static_cast<int>( static_cast<float>(SizeX - 1) / AspectRatio) + 1;
	}
}

