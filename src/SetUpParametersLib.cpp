/************************************************************************************************\
* 2020 Pedro Cosme , João Santos, Ivan Figueiredom, João Rebelo, Diogo Simões                    *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/

#include "includes/TethysBaseLib.h"
#include "includes/SetUpParametersLib.h"

SetUpParameters::SetUpParameters(){
	SizeX=101;
	SizeY=101;
	SoundVelocity = 30.0f;
	FermiVelocity = 10.0f;
	CollisionFrequency = 0.01f;
	ShearViscosity = 0.0f;
	CyclotronFrequency = 0.0f;
	OddViscosity=0.0f;
	ThermalDiffusivity=0.0f;
	SaveMode = 1;
	this->DefineGeometry();
}

SetUpParameters::SetUpParameters(float sound, float fermi, float coll, float visco, float odd_visco, float cyclo,
								 float thermal_diff, int mode, float aspect) {
	SizeX=101;
	SizeY=101;
	SoundVelocity = sound;
	FermiVelocity = fermi;
	CollisionFrequency = coll;
	ShearViscosity = visco;
	CyclotronFrequency = cyclo;
	SaveMode = mode;
	AspectRatio = aspect;
	OddViscosity=odd_visco;
	ThermalDiffusivity=thermal_diff;
	ParametersChecking();
	DefineGeometry();
}

SetUpParameters::SetUpParameters(int argc, char ** argv) {
	SizeX=101;
	SizeY=101;

	if(argc==2){
		ReadIniFile(argv[1]);
	}

	else{
		if (argc == 9 || argc == 10) {
			SoundVelocity = strtof(argv[1], nullptr);
			FermiVelocity = strtof(argv[2], nullptr);
			CollisionFrequency = strtof(argv[3], nullptr);
			ShearViscosity = strtof(argv[4], nullptr);
			OddViscosity = strtof(argv[5], nullptr);
			ThermalDiffusivity = strtof(argv[6], nullptr);
			CyclotronFrequency = strtof(argv[7], nullptr);
			SaveMode = (int) strtol(argv[8], nullptr, 10);    // full data or light save option
			if (argc == 10) {
				AspectRatio = strtof(argv[9], nullptr);
			}
		}
		else {
			this->PromptParameters();
		}
	}
	ParametersChecking();
	DefineGeometry();
}

void SetUpParameters::ParametersChecking() const{
	std::string msg;
	if(SoundVelocity<=0.0f){
		msg = "ERROR: Unphysical Sound Velocity";
		cerr << msg <<"\nExiting"<< endl;
		exit(EXIT_FAILURE);
	}
	if(FermiVelocity<0.0f){
		msg = "ERROR: Unphysical Fermi Velocity";
		cerr << msg <<"\nExiting"<< endl;
		exit(EXIT_FAILURE);
	}
	if(ShearViscosity<0.0f){
		msg = "ERROR: Unphysical Shear Viscosity";
		cerr << msg <<"\nExiting"<< endl;
		exit(EXIT_FAILURE);
	}
	if(OddViscosity<0.0f){
		msg = "ERROR: Unphysical Odd Viscosity";
		cerr << msg <<"\nExiting"<< endl;
		exit(EXIT_FAILURE);
	}
	if(CollisionFrequency<0.0f){
		msg = "ERROR: Unphysical Collision Frequency";
		cerr << msg <<"\nExiting"<< endl;
		exit(EXIT_FAILURE);
	}
	if(ThermalDiffusivity<0.0f){
		msg = "ERROR: Thermal Diffusivity";
		cerr << msg <<"\nExiting"<< endl;
		exit(EXIT_FAILURE);
	}
	if(CyclotronFrequency<0.0f){
		msg = "ERROR: Unphysical Cyclotron Frequency";
		cerr << msg <<"\nExiting"<< endl;
		exit(EXIT_FAILURE);
	}
	if(AspectRatio<0.0f){
		msg = "ERROR: Negative Aspect Ratio";
		cerr << msg <<"\nExiting"<< endl;
		exit(EXIT_FAILURE);
	}
	if( SaveMode != 0 && SaveMode != 1  ) {
		msg = "ERROR: Unknown save mode option";
		cerr << msg <<"\nExiting"<< endl;
		exit(EXIT_FAILURE);
	}

}

void SetUpParameters::ParametersFromHdf5File(const string& hdf5name){
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
	Attribute attr_vis(grp_dat->openAttribute("Kinematic shear viscosity"));
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
	this->DefineGeometry();
}

void SetUpParameters::DefineGeometry() {
	int sampling=201;
	if(AspectRatio>1.0f){
		Length=1.0f*AspectRatio;
		Width=1.0f;
		SizeY=sampling;
		SizeX= static_cast<int>( static_cast<float>(SizeY-1)*AspectRatio)+1;
	}
	if(AspectRatio==1.0f){
		Length=1.0f;
		Width=1.0f;
		SizeX=sampling;
		SizeY=sampling;
	}
	if(AspectRatio<1.0f){
		Length=1.0f;
		Width=1.0f/AspectRatio;
		SizeX=sampling;
		SizeY= static_cast<int>( static_cast<float>(SizeX - 1) / AspectRatio) + 1;
	}
}


void SetUpParameters::PrintParameters() const{
	cout << endl << "Sound velocity: " << SoundVelocity << endl;
	cout << "Fermi velocity: " << FermiVelocity << endl;
	cout << "Shear viscosity: " << ShearViscosity << endl;
	cout << "Odd viscosity: " << OddViscosity << endl;
	cout << "Collision frequency: " << CollisionFrequency << endl;
	cout << "Cyclotron frequency: " << CyclotronFrequency << endl;
	cout << "Thermal Diffusivity: " << ThermalDiffusivity << endl;
	cout << "Aspect ratio: " << AspectRatio << endl;
	cout << "Save mode: " << SaveMode << endl;
}

void SetUpParameters::ReadIniFile(char *  file_name) {
	string line;
	ifstream file(file_name);

	vector<string> veclines;

	while(getline(file, line)){
		veclines.push_back(line);
	}
	file.close();

	for(auto & vecline : veclines){
		line = vecline;
		vector<char> num_part, txt_part;
		if(line[0]!=';' && line[0]!='[') {
			for (char &i : line) {
				if (isdigit(i) || int(i) == 46)
					num_part.push_back(i);
				else {
					if ((int(i) >= 65 && int(i) <= 90) || (int(i) >= 97 && int(i) <= 122))
						txt_part.push_back(i);
				}
			}
		}
		string num(num_part.begin(), num_part.end());
		string txt(txt_part.begin(), txt_part.end());

		if(txt == "sound") SoundVelocity = stof(num);
		if(txt == "fermi") FermiVelocity = stof(num);
		if(txt == "shear") ShearViscosity = stof(num);
		if(txt == "odd") OddViscosity = stof(num);
		if(txt == "col") CollisionFrequency = stof(num);
		if(txt == "cycl") CyclotronFrequency = stof(num);
		if(txt == "therm") ThermalDiffusivity = stof(num);
		if(txt == "aspect") AspectRatio = stof(num);
		if(txt == "time") SimulationTime = stof(num);
		if(txt == "save") SaveMode = stoi(num);
	}
}

void SetUpParameters::PromptParameters() {
	cout << "Define S value: ";
	cin >> SoundVelocity;
	cout << "Define vF value: ";
	cin >> FermiVelocity;
	cout << "Define kinetic shear viscosity: ";
	cin >> ShearViscosity;
	cout << "Define kinetic odd viscosity: ";
	cin >> OddViscosity;
	cout << "Define collision frequency: ";
	cin >> CollisionFrequency;
	cout << "Define cyclotron frequency: ";
	cin >> CyclotronFrequency;
	cout << "Define thermal diffusivity: ";
	cin >> ThermalDiffusivity;
	cout << "Define the aspect ratio x:y ";
	cin >> AspectRatio;
	cout << "Define data_save_mode value (0-> light save | 1-> full data): ";
	cin >> SaveMode;
}


void SetUpParametersCNP::PromptParameters() {
	cout << "Define S value: ";
	cin >> SoundVelocity;
	cout << "Define vF value: ";
	cin >> FermiVelocity;
	cout << "Define vT value (Dirac): ";
	cin >> ThermalVelocity;
	cout << "Define A source term value (Dirac): ";
	cin >> Diffusive_sourceterm;
	cout << "Define B source term value (Dirac): ";
	cin >> Creation_sourceterm;
	cout << "Define kinetic shear viscosity: ";
	cin >> ShearViscosity;
	cout << "Define kinetic odd viscosity: ";
	cin >> OddViscosity;
	cout << "Define collision frequency: ";
	cin >> CollisionFrequency;
	cout << "Define cyclotron frequency: ";
	cin >> CyclotronFrequency;
	cout << "Define thermal diffusivity: ";
	cin >> ThermalDiffusivity;
	cout << "Define the aspect ratio x:y ";
	cin >> AspectRatio;
	cout << "Define data_save_mode value (0-> light save | 1-> full data): ";
	cin >> SaveMode;
}


void SetUpParametersCNP::ReadIniFile(char *  file_name) {
	string line;
	ifstream file(file_name);

	vector<string> veclines;

	while(getline(file, line)){
		veclines.push_back(line);
	}
	file.close();

	for(auto & vecline : veclines){
		line = vecline;
		vector<char> num_part, txt_part;
		if(line[0]!=';' && line[0]!='[') {
			for (char &i : line) {
				if (isdigit(i) || int(i) == 46)
					num_part.push_back(i);
				else {
					if ((int(i) >= 65 && int(i) <= 90) || (int(i) >= 97 && int(i) <= 122))
						txt_part.push_back(i);
				}
			}
		}
		string num(num_part.begin(), num_part.end());
		string txt(txt_part.begin(), txt_part.end());

		if(txt == "vth,") ThermalVelocity = stof(num);
		if(txt == "sourceDiff") Diffusive_sourceterm = stof(num);
		if(txt == "sourceCrea") Creation_sourceterm = stof(num);

		if(txt == "sound") SoundVelocity = stof(num);
		if(txt == "fermi") FermiVelocity = stof(num);
		if(txt == "shear") ShearViscosity = stof(num);
		if(txt == "odd") OddViscosity = stof(num);
		if(txt == "col") CollisionFrequency = stof(num);
		if(txt == "cycl") CyclotronFrequency = stof(num);
		if(txt == "therm") ThermalDiffusivity = stof(num);
		if(txt == "aspect") AspectRatio = stof(num);
		if(txt == "time") SimulationTime = stof(num);
		if(txt == "save") SaveMode = stoi(num);
	}
}


SetUpParametersCNP::SetUpParametersCNP(int argc, char ** argv) {
	SizeX=101;
	SizeY=101;

	if(argc==2){
		this->ReadIniFile(argv[1]);
	}

	else{
		if (argc == 9 || argc == 10) {
			SoundVelocity = strtof(argv[1], nullptr);
			FermiVelocity = strtof(argv[2], nullptr);
			CollisionFrequency = strtof(argv[3], nullptr);
			ShearViscosity = strtof(argv[4], nullptr);
			OddViscosity = strtof(argv[5], nullptr);
			ThermalDiffusivity = strtof(argv[6], nullptr);
			CyclotronFrequency = strtof(argv[7], nullptr);
			SaveMode = (int) strtol(argv[8], nullptr, 10);    // full data or light save option
			if (argc == 10) {
				AspectRatio = strtof(argv[9], nullptr);
			}
		}
		else {
			this->PromptParameters();
		}
	}
	this->ParametersChecking();
	DefineGeometry();
}


SetUpParametersCNP::SetUpParametersCNP(){
	SizeX=101;
	SizeY=101;
	SoundVelocity = 30.0f;
	FermiVelocity = 10.0f;
	CollisionFrequency = 0.01f;
	ShearViscosity = 0.0f;
	CyclotronFrequency = 0.0f;
	OddViscosity=0.0f;
	ThermalDiffusivity=0.0f;
	SaveMode = 1;
	this->DefineGeometry();
}

/*
void SetUpParameters1D::PromptParameters() {
	cout << "Define S value: ";
	cin >> SoundVelocity;
	cout << "Define vF value: ";
	cin >> FermiVelocity;
	cout << "Define kinetic shear viscosity: ";
	cin >> ShearViscosity;
	cout << "Define kinetic odd viscosity: ";
	cin >> CollisionFrequency;
	cout << "Define thermal diffusivity: ";
	cin >> ThermalDiffusivity;
	cout << "Define data_save_mode value (0-> light save | 1-> full data): ";
	cin >> SaveMode;
}

void SetUpParameters1D::ReadIniFile(char *  file_name) {
	string line;
	ifstream file(file_name);

	vector<string> veclines;

	while(getline(file, line)){
		veclines.push_back(line);
	}
	file.close();

	for(auto & vecline : veclines){
		line = vecline;
		vector<char> num_part, txt_part;
		if(line[0]!=';' && line[0]!='[') {
			for (char &i : line) {
				if (isdigit(i) || int(i) == 46)
					num_part.push_back(i);
				else {
					if ((int(i) >= 65 && int(i) <= 90) || (int(i) >= 97 && int(i) <= 122))
						txt_part.push_back(i);
				}
			}
		}
		string num(num_part.begin(), num_part.end());
		string txt(txt_part.begin(), txt_part.end());

		if(txt == "sound") SoundVelocity = stof(num);
		if(txt == "fermi") FermiVelocity = stof(num);
		if(txt == "shear") ShearViscosity = stof(num);
		if(txt == "col") CollisionFrequency = stof(num);
		if(txt == "therm") ThermalDiffusivity = stof(num);
		if(txt == "time") SimulationTime = stof(num);
		if(txt == "save") SaveMode = stoi(num);
	}
}
 */