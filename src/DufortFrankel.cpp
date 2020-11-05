//
// Created by pcosme on 05/11/2020.
//

#include "Tethys2DLib.h"
#include "BoundaryLib.h"
#include "ElectricLib.h"

#ifndef MAT_PI
#	define MAT_PI 3.14159265358979323846
#endif


using namespace std;



int main(int argc, char **argv){
	SetUpParameters parametros(30.0f,10.0f,0.0f,5.0f,0.0f,1,1.0);
	parametros.DefineGeometry();

	GrapheneFluid2D grafeno(parametros);
	DirichletBoundaryCondition boundary_condition;

	float t=0.0;
	float dt;		// time step

	/*......CFL routine to determine dt...............................*/
	grafeno.CflCondition();
	dt=grafeno.GetDt();
	grafeno.SetTmax(.2f);

	/*.........Output files and streams...............................*/
	grafeno.CreateFluidFile();
	grafeno.CreateHdf5File();
	if(parametros.SaveMode){
		grafeno.SaveSound();
	}
	//grafeno.InitialCondRand();
	grafeno.InitialCondTest();

	/*................................................................*/
	grafeno.SetSnapshotStep(5);
	while (t <= grafeno.GetTmax() ){

		t += dt;
		GrapheneFluid2D::TimeStepCounter++;
		cout<<t<<"\t"<<GrapheneFluid2D::TimeStepCounter<<endl;


		//grafeno.ParabolicOperatorFtcs();
		grafeno.ParabolicOperatorDuFortFrankel();

		grafeno.TimeUpdate();

		boundary_condition.YFree(grafeno);
		boundary_condition.XFree(grafeno);



		//Record full hdf5 data
		if (parametros.SaveMode  && grafeno.Snapshot()) {
			grafeno.SaveSnapShot();
		}
		grafeno.WriteFluidFile(t);
	}
	//Record atributes on hdf5 file
	if(parametros.SaveMode) {
		grafeno.WriteAttributes();
	}
	grafeno.CloseHdf5File();

	return 0;
}