#include "includes/Tethys2DLib.h"
#include "includes/BoundaryLib.h"
#include "includes/ElectricLib.h"

#ifndef MAT_PI
#	define MAT_PI 3.14159265358979323846
#endif


using namespace std;


int main(int argc, char **argv){
	SetUpParameters parameters(argc, argv);
	parameters.DefineGeometry();

	float t=0.0;
	float dt;		// time step

	GrapheneFluid2D graph(parameters);

	/*......CFL routine to determine dt...............................*/
	graph.CflCondition();
	dt=graph.GetDt();
	/*................................................................*/
	
	/*.........Fixed or variable vel_snd value........................*/
	graph.SetSound();
	//graph.SetSimulationTime();
	graph.SetTmax(3.0f);
	/*................................................................*/

	/*.........Output files and streams...............................*/
	graph.CreateFluidFile();
	graph.CreateHdf5File();
	if(parameters.SaveMode){
		graph.SaveSound();
	}
	/*................................................................*/

	GrapheneFluid2D::BannerDisplay();
	graph.WelcomeScreen();

	/*...............Initialization...................................*/
	graph.InitialCondRand();
	/*................................................................*/

	/*................Setting.the.lateral.boundaries..................*/
	DyakonovShurBoundaryCondition::SetSlope(0.0f);
	DyakonovShurBoundaryCondition::SetBottomEdge(graph);
	DyakonovShurBoundaryCondition::SetTopEdge(graph);
	/*................................................................*/


	cout << "\033[1;7;5;33m Program Running \033[0m"<<endl;
	while (t <= graph.GetTmax() ){

		t += dt;
		GrapheneFluid2D::TimeStepCounter++;

		graph.Richtmyer();
		DyakonovShurBoundaryCondition::DyakonovShurBc(graph);
		DyakonovShurBoundaryCondition::YFree(graph);

		if(graph.GetKinVis()!=0.0f ) {
			graph.ParabolicOperatorWeightedExplicit19();
			DyakonovShurBoundaryCondition::DyakonovShurBc(graph);
			DyakonovShurBoundaryCondition::YFree(graph);
		}

		//Record full hdf5 data
		if (parameters.SaveMode  && graph.Snapshot()) {
			graph.SaveSnapShot();
		}
		if(static_cast<int>(fmod(t/dt,2.0f))){
			graph.WriteFluidFile(t);
		}
	}
	//Record atributes on hdf5 file
	if(parameters.SaveMode) {
		graph.WriteAttributes();
	}
	graph.CloseHdf5File();


	cout << "\033[1A\033[2K\033[1;32mDONE!\033[0m\n";
	cout << "═══════════════════════════════════════════════════════════════════════════" <<endl;
	return 0;
}