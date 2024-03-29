/************************************************************************************************\
* 2020 Pedro Cosme , João Santos, Ivan Figueiredom, João Rebelo, Diogo Simões                    *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/

#include "includes/ElectricLib.h"
#include "includes/GrapheneFluid2DLib.h"
#include "includes/GrapheneFluid1DLib.h"

#ifndef MAT_PI
#    define MAT_PI 3.14159265358979323846
#endif

#ifndef MAT_EULER
#    define MAT_EULER 2.71828182845905
#endif

using namespace H5;
using namespace std;


void ElectroAnalysis::CreateElectroFile(const GrapheneFluid1D& graphene){
	std::string infix = graphene.GetInfix();
	std::string electrofile;
	electrofile = "electro_1D_" + infix + ".dat" ;
	data_electro.open (electrofile);
	data_electro << scientific;
}

void ElectroAnalysis::CreateElectroFile(const string &infix_string){
	std::string electrofile;
	electrofile = "electro_2D_" + infix_string + ".dat" ;
	data_electro.open (electrofile);
	data_electro << scientific;
}

void ElectroAnalysis::WriteElectroFile(float t,const GrapheneFluid1D& graphene){
	float q_net = ElectroAnalysis::NetCharge(graphene);
	float i_avg = ElectroAnalysis::AverageCurrent(graphene);
	float p_ohm = ElectroAnalysis::OhmPower(graphene);
	float dipole_var=ElectroAnalysis::ElectricDipoleVariation(graphene);
	float dipole=ElectroAnalysis::ElectricDipole(graphene);
	data_electro << t << "\t" << q_net << "\t" << i_avg << "\t" << q_net * q_net * 0.5 << "\t" << p_ohm << "\t" << dipole << "\t" << dipole_var << "\n";
}

void ElectroAnalysis::ComputeElectroBase(float t, const GrapheneFluid2D& graphene){
	TmpArr.push_back(t);
	NetQ.push_back(ElectroAnalysis::NetCharge(graphene));
	DipX.push_back(ElectroAnalysis::ElectricDipoleX(graphene));
	DipY.push_back(ElectroAnalysis::ElectricDipoleY(graphene));
	QuadXX.push_back(ElectroAnalysis::ElectricQuadrupoleXX(graphene));
	QuadXY.push_back(ElectroAnalysis::ElectricQuadrupoleXY(graphene));
	QuadYY.push_back(ElectroAnalysis::ElectricQuadrupoleYY(graphene));
	DipMagZ.push_back(ElectroAnalysis::MagneticDipoleZ(graphene));
	CurD.push_back(ElectroAnalysis::DrainCurrent(graphene));
	CurS.push_back(ElectroAnalysis::SourceCurrent(graphene));
	AvgCurDS.push_back(ElectroAnalysis::AverageDirectCurrent(graphene));
	VoltDS.push_back(ElectroAnalysis::DrainToSourceVoltage(graphene));
	AvgCurHall.push_back(ElectroAnalysis::AverageHallCurrent(graphene));
	PowOhm.push_back(ElectroAnalysis::OhmPower(graphene));
}


void ElectroAnalysis::ComputeElectroDerived() {
	float dt;
	dt=TmpArr.back()/static_cast<float>(DipX.size());
	EngCap.resize(DipX.size());
	PowCap.resize(DipX.size());
	transform(NetQ.begin(), NetQ.end(), EngCap.begin(), [](const float &c){ return 0.5f*c*c; });
	Convolve_Gauss(1,5,1.0,EngCap.data(),PowCap.data(),EngCap.size());
	transform(PowCap.begin(), PowCap.end(), PowCap.begin(), [dt](const float &c){ return c/dt; });
	DipVarX.resize(DipX.size());
	DipVarVarX.resize(DipX.size());
	DipVarY.resize(DipY.size());
	DipVarVarY.resize(DipY.size());
	QuadVarXX.resize(QuadXX.size());
	QuadVarVarXX.resize(QuadXX.size());
	QuadVarVarVarXX.resize(QuadXX.size());
	QuadVarXY.resize(QuadXY.size());
	QuadVarVarXY.resize(QuadXY.size());
	QuadVarVarVarXY.resize(QuadXY.size());
	QuadVarYY.resize(QuadYY.size());
	QuadVarVarYY.resize(QuadYY.size());
	QuadVarVarVarYY.resize(QuadYY.size());
	DipMagVarZ.resize(DipMagZ.size());
	DipMagVarVarZ.resize(DipMagZ.size());
	Convolve_Gauss(1,5,1.0,DipX.data(),DipVarX.data(),DipX.size());
	Convolve_Gauss(1,5,1.0,DipVarX.data(),DipVarVarX.data(),DipX.size());
	Convolve_Gauss(1,5,1.0,DipY.data(),DipVarY.data(),DipY.size());
	Convolve_Gauss(1,5,1.0,DipVarY.data(),DipVarVarY.data(),DipY.size());
	Convolve_Gauss(1,5,1.0,QuadXX.data(),QuadVarXX.data(),QuadXX.size());
	Convolve_Gauss(1,5,1.0,QuadVarXX.data(),QuadVarVarXX.data(),QuadXX.size());
	Convolve_Gauss(1,5,1.0,QuadVarVarXX.data(),QuadVarVarVarXX.data(),QuadXX.size());
	Convolve_Gauss(1,5,1.0,QuadXY.data(),QuadVarXY.data(),QuadXY.size());
	Convolve_Gauss(1,5,1.0,QuadVarXY.data(),QuadVarVarXY.data(),QuadXY.size());
	Convolve_Gauss(1,5,1.0,QuadVarVarXY.data(),QuadVarVarVarXY.data(),QuadXY.size());
	Convolve_Gauss(1,5,1.0,QuadYY.data(),QuadVarYY.data(),QuadYY.size());
	Convolve_Gauss(1,5,1.0,QuadVarYY.data(),QuadVarVarYY.data(),QuadYY.size());
	Convolve_Gauss(1,5,1.0,QuadVarVarYY.data(),QuadVarVarVarYY.data(),QuadYY.size());
	Convolve_Gauss(1,5,1.0,DipMagZ.data(),DipMagVarZ.data(),DipMagZ.size());
	Convolve_Gauss(1,5,1.0,DipMagVarZ.data(),DipMagVarVarZ.data(),DipMagZ.size());
	transform(DipVarX.begin(), DipVarX.end(), DipVarX.begin(), [dt](const float &c){ return c/dt; });
	transform(DipVarVarX.begin(), DipVarVarX.end(), DipVarVarX.begin(), [dt](const float &c){ return c/(dt*dt); });
	transform(DipVarY.begin(), DipVarY.end(), DipVarY.begin(), [dt](const float &c){ return c/dt; });
	transform(DipVarVarY.begin(), DipVarVarY.end(), DipVarVarY.begin(), [dt](const float &c){ return c/(dt*dt); });
	transform(QuadVarXX.begin(), QuadVarXX.end(), QuadVarXX.begin(), [dt](const float &c){ return c/dt; });
	transform(QuadVarVarXX.begin(), QuadVarVarXX.end(), QuadVarVarXX.begin(), [dt](const float &c){ return c/(dt*dt); });
	transform(QuadVarVarVarXX.begin(), QuadVarVarVarXX.end(), QuadVarVarVarXX.begin(), [dt](const float &c){ return c/(dt*dt*dt); });
	transform(QuadVarXY.begin(), QuadVarXY.end(), QuadVarXY.begin(), [dt](const float &c){ return c/dt; });
	transform(QuadVarVarXY.begin(), QuadVarVarXY.end(), QuadVarVarXY.begin(), [dt](const float &c){ return c/(dt*dt); });
	transform(QuadVarVarVarXY.begin(), QuadVarVarVarXY.end(), QuadVarVarVarXY.begin(), [dt](const float &c){ return c/(dt*dt*dt); });
	transform(QuadVarYY.begin(), QuadVarYY.end(), QuadVarYY.begin(), [dt](const float &c){ return c/dt; });
	transform(QuadVarVarYY.begin(), QuadVarVarYY.end(), QuadVarVarYY.begin(), [dt](const float &c){ return c/(dt*dt); });
	transform(QuadVarVarVarYY.begin(), QuadVarVarVarYY.end(), QuadVarVarVarYY.begin(), [dt](const float &c){ return c/(dt*dt*dt); });
	transform(DipMagVarZ.begin(), DipMagVarZ.end(), DipMagVarZ.begin(), [dt](const float &c){ return c/dt; });
	transform(DipMagVarVarZ.begin(), DipMagVarVarZ.end(), DipMagVarVarZ.begin(), [dt](const float &c){ return c/(dt*dt); });
}


void ElectroAnalysis::WriteElectroFile() {
		if(NetQ.empty() || DipVarY.empty()){
			cerr << "Nothing to save on output file" <<"\nExiting"<< endl;
			exit(EXIT_FAILURE);
		}
		for(size_t i = 0; i < NetQ.size(); ++i) {
			data_electro << TmpArr[i] << "\t"
			             << NetQ[i] << "\t"
			             << AvgCurDS[i] << "\t"
			             << AvgCurHall[i] << "\t"
			             << VoltDS[i] << "\t"
			             << CurD[i] << "\t"
			             << CurS[i] << "\t"
			             << PowOhm[i] << "\t"
			             << PowCap[i] << "\t"
			             << DipX[i] << "\t"
			             << DipVarX[i] << "\t"
			             << DipVarVarX[i] << "\t"
			             << DipY[i] << "\t"
			             << DipVarY[i] << "\t"
			             << DipVarVarY[i] << "\t"
			             << QuadXX[i] << "\t"
			             << QuadVarXX[i] << "\t"
			             << QuadVarVarXX[i] << "\t"
			             << QuadVarVarVarXX[i] << "\t"
			             << QuadXY[i] << "\t"
			             << QuadVarXY[i] << "\t"
			             << QuadVarVarXY[i] << "\t"
			             << QuadVarVarVarXY[i] << "\t"
			             << QuadYY[i] << "\t"
			             << QuadVarYY[i] << "\t"
			             << QuadVarVarYY[i] << "\t"
			             << QuadVarVarVarYY[i] << "\t"
			             << DipMagZ[i] << "\t"
			             << DipMagVarZ[i] << "\t"
			             << DipMagVarVarZ[i] << "\n";
		}
}

float ElectroAnalysis::NetCharge(const GrapheneFluid1D& graphene){
	return Integral_1_D(graphene.SizeX(), graphene.GetDx(), graphene.Den);
}
float ElectroAnalysis::NetCharge(const GrapheneFluid2D& graphene){
	return Integral_2_D(graphene.SizeX(),graphene.SizeY(), graphene.GetDx(),graphene.GetDy(), graphene.Den);
}

float ElectroAnalysis::AverageCurrent(const GrapheneFluid1D& graphene){
	return Integral_1_D(graphene.SizeX(), graphene.GetDx(), graphene.Den);
}

float ElectroAnalysis::AverageDirectCurrent(const GrapheneFluid2D& graphene){
	return Integral_2_D(graphene.SizeX(),graphene.SizeY(), graphene.GetDx(),graphene.GetDy(), graphene.CurX);
}

float ElectroAnalysis::AverageHallCurrent(const GrapheneFluid2D& graphene){
	return Integral_2_D(graphene.SizeX(),graphene.SizeY(), graphene.GetDx(),graphene.GetDy(), graphene.CurY);
}


float ElectroAnalysis::ElectricDipoleVariation(const GrapheneFluid1D& graphene){
	return Integral_1_D(graphene.SizeX(), graphene.GetDx(), graphene.Den);
}

float ElectroAnalysis::ElectricDipole(const GrapheneFluid1D& graphene){
	float p=0.0;
	float dx=graphene.GetDx();
	for(int j=1;j<graphene.SizeX()/2;j++){
		p += dx*static_cast<float>(2*j-2)*graphene.Den[2 * j - 2] + 4 * dx * static_cast<float>(2 * j - 1) * graphene.Den[2 * j - 1] + dx * static_cast<float>(2 * j) * graphene.Den[2 * j];
	}
	p = p*graphene.GetDx()/3.0f;
	return p;
}

float ElectroAnalysis::OhmPower(const GrapheneFluid1D& graphene){
	float itg=0.0;
	for(int j=1;j<graphene.SizeX()/2;j++){
		itg += graphene.Den[2 * j - 2] * graphene.Den[2 * j - 2]/sqrt(graphene.Den[2 * j - 2]) + 4 * graphene.Den[2 * j - 1] * graphene.Den[2 * j - 1] /sqrt(graphene.Den[2 * j - 1]) + graphene.Den[2 * j] * graphene.Den[2 * j]/sqrt(graphene.Den[2 * j]);
	}
	itg = itg*graphene.GetDx()/3.0f;
	return itg;
}

float ElectroAnalysis::OhmPower(const GrapheneFluid2D& graphene){
	int size = graphene.SizeX()*graphene.SizeY();
	float square_current_density[size];
	for(int c=0;c<size;c++){
		float jx,jy;
		jx=graphene.CurX[c];
		jy=graphene.CurY[c];
		square_current_density[c] = (jx*jx+jy*jy)/sqrt(graphene.Den[c]);
	}
	return Integral_2_D(graphene.SizeX(), graphene.SizeY(), graphene.GetDx(), graphene.GetDy(), square_current_density);
}

float ElectroAnalysis::ElectricDipoleX(const GrapheneFluid2D &graphene) {
	int size = graphene.SizeX()*graphene.SizeY();
	float vector[size];
	for(int c=0;c<size;c++){
		div_t divresult;
		divresult = div (c,graphene.SizeX());
		auto i=static_cast<float>(divresult.rem);
		float rx;
		rx = i*graphene.GetDx()-0.5f*graphene.GetLengthX();
		vector[c] = rx*graphene.Den[c];
	}
	return Integral_2_D(graphene.SizeX(), graphene.SizeY(), graphene.GetDx(), graphene.GetDy(), vector);
}
float ElectroAnalysis::ElectricDipoleY(const GrapheneFluid2D &graphene) {
	int size = graphene.SizeX()*graphene.SizeY();
	float vector[size];
	for(int c=0;c<size;c++){
		div_t divresult;
		divresult = div (c,graphene.SizeX());
		auto j=static_cast<float>(divresult.quot);
		float ry;
		ry = j*graphene.GetDy()-0.5f*graphene.GetLengthY();
		vector[c] = ry*graphene.Den[c];
	}
	return Integral_2_D(graphene.SizeX(), graphene.SizeY(), graphene.GetDx(), graphene.GetDy(), vector);
}
float ElectroAnalysis::ElectricDipoleVariationX(const GrapheneFluid2D &graphene) {
	return Integral_2_D(graphene.SizeX(), graphene.SizeY(), graphene.GetDx(), graphene.GetDy(), graphene.CurX );
}
float ElectroAnalysis::ElectricDipoleVariationY(const GrapheneFluid2D &graphene) {
	return Integral_2_D(graphene.SizeX(), graphene.SizeY(), graphene.GetDx(), graphene.GetDy(), graphene.CurY );
}

float ElectroAnalysis::ElectricQuadrupoleXX(const GrapheneFluid2D &graphene){
	int size = graphene.SizeX()*graphene.SizeY();
	float vector[size];
	for(int c=0;c<size;c++){
		div_t divresult;
		divresult = div (c,graphene.SizeX());
		auto i=static_cast<float>(divresult.rem);
		auto j=static_cast<float>(divresult.quot);
		float rx;
		rx = i*graphene.GetDx()-0.5f*graphene.GetLengthX();
		float ry;
		ry = j*graphene.GetDy()-0.5f*graphene.GetLengthY();
		vector[c] = (2*rx*rx-ry*ry)*graphene.Den[c];
	}
	return Integral_2_D(graphene.SizeX(), graphene.SizeY(), graphene.GetDx(), graphene.GetDy(), vector);
}

float ElectroAnalysis::ElectricQuadrupoleXY(const GrapheneFluid2D &graphene){
	int size = graphene.SizeX()*graphene.SizeY();
	float vector[size];
	for(int c=0;c<size;c++){
		div_t divresult;
		divresult = div (c,graphene.SizeX());
		auto i=static_cast<float>(divresult.rem);
		auto j=static_cast<float>(divresult.quot);
		float rx;
		rx = i*graphene.GetDx()-0.5f*graphene.GetLengthX();
		float ry;
		ry = j*graphene.GetDy()-0.5f*graphene.GetLengthY();
		vector[c] = (2*ry*ry-rx*rx)*graphene.Den[c];
	}
	return Integral_2_D(graphene.SizeX(), graphene.SizeY(), graphene.GetDx(), graphene.GetDy(), vector);
}

float ElectroAnalysis::ElectricQuadrupoleYY(const GrapheneFluid2D &graphene){
	int size = graphene.SizeX()*graphene.SizeY();
	float vector[size];
	for(int c=0;c<size;c++){
		div_t divresult;
		divresult = div (c,graphene.SizeX());
		auto i=static_cast<float>(divresult.rem);
		auto j=static_cast<float>(divresult.quot);
		float rx;
		rx = i*graphene.GetDx()-0.5f*graphene.GetLengthX();
		float ry;
		ry = j*graphene.GetDy()-0.5f*graphene.GetLengthY();
		vector[c] = 3*rx*ry*graphene.Den[c];
	}
	return Integral_2_D(graphene.SizeX(), graphene.SizeY(), graphene.GetDx(), graphene.GetDy(), vector);
}
//\int_0^{W/L}\int_0^1 ( (x-w/2)*j_y - (y-1/2)*j_x ) \,dxdy
float ElectroAnalysis::MagneticDipoleZ(const GrapheneFluid2D &graphene){
	int size = graphene.SizeX()*graphene.SizeY();
	float vector[size];
	for(int c=0;c<size;c++){
		div_t divresult;
		divresult = div (c,graphene.SizeX());
		auto i=static_cast<float>(divresult.rem);
		auto j=static_cast<float>(divresult.quot);
		float rx;
		rx = i*graphene.GetDx()-0.5f*graphene.GetLengthX();
		float ry;
		ry = j*graphene.GetDy()-0.5f*graphene.GetLengthY();
		vector[c] = rx*graphene.CurY[c]-ry*graphene.CurX[c];
	}
	return Integral_2_D(graphene.SizeX(), graphene.SizeY(), graphene.GetDx(), graphene.GetDy(), vector);
}

void ElectroAnalysis::CloseElectroFile() {
	data_electro.close();
}

float ElectroAnalysis::DrainToSourceVoltage(const GrapheneFluid2D& graphene) {
	return Integral_2_D(graphene.SizeX(), graphene.SizeY(), graphene.GetDx(), graphene.GetDy(), graphene.VelX );
}

float ElectroAnalysis::SourceCurrent(const GrapheneFluid2D& graphene) {
	int size1d = graphene.SizeY();
	float vector[size1d];
	for(int j=0;j<size1d;j++){
		int pos;
		pos=j*graphene.SizeX();
		vector[j] = graphene.CurX[pos];
	}
	return Integral_1_D(graphene.SizeY(), graphene.GetDy(), vector);
}

float ElectroAnalysis::DrainCurrent(const GrapheneFluid2D& graphene) {
	int size1d = graphene.SizeY();
	float vector[size1d];
	for(int j=0;j<size1d;j++){
		int pos;
		pos=(j+1)*graphene.SizeX()-1;
		vector[j] = graphene.CurX[pos];
	}
	return Integral_1_D(graphene.SizeY(), graphene.GetDy(), vector);
}

void ElectroAnalysis::BannerDisplay(const GrapheneFluid2D &graphene) {
	GrapheneFluid2D::BannerDisplay();
	cout<<"┌─────────────────────────────────────────────────────────────────────────┐\n";
	cout<<"│                \033[3mComputation of the Electronic Quantities\033[0m                 │\n";
	cout<<"└─────────────────────────────────────────────────────────────────────────┘\n";

	cout<<"Detected parameters: \n";
	cout << "Fermi velocity\t\033[1mvF\t" << graphene.GetVelFer() << " v\342\202\200\033[0m\n";
	cout << "Phase velocity\t\033[1mS'\t" << graphene.PhaseVel() << " v\342\202\200\033[0m\n";
	cout << "Viscosity \t\033[1m\316\267\t" << graphene.GetKinVis() << "\033[0m\n";
	if (graphene.GetKinVis() != 0.0) {
		cout << "Reynolds n. \t\033[1mRe\t" << 1.0 / graphene.GetKinVis() << "\033[0m\n";
	}
	cout << "Collision rate \t\033[1m\316\275\t" << graphene.GetColFreq()<< " v\342\202\200/L\033[0m\n";
	cout << "Cyclotron frequency \t\033[1m\317\211c\t" << graphene.GetCycFreq() << " v\342\202\200/L\n\033[0m\n";
}

