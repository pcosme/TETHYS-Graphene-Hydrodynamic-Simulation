#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <string>
#include <ctime>
#include <numeric>
//#include <fftw3.h>

#include "Tethys1DLib.h"

using namespace std;




int main(int argc, char **argv){
	cout<<"\n" ;
	cout<<"╔═════════════════════════════════════════════════════════════════════════╗\n";
	cout<<"║ \033[1m Electronic analysis for TETHYS                                        \033[0m ║\n";
	cout<<"╚═════════════════════════════════════════════════════════════════════════╝\n";       	

	int N = atoi(argv[1]);	
	
	float S = atof(argv[3]);	
	
	float *Time;			
	Time =(float*) calloc (N,sizeof(float));
	
	float *in_net_Q;			
	in_net_Q =(float*) calloc (N,sizeof(float));


	float * in_avg_I;
	in_avg_I=(float*) calloc (N,sizeof(float));

	float *in_sto_E;			
	in_sto_E =(float*) calloc (N,sizeof(float));

	float *in_P_Ohm;			
	in_P_Ohm =(float*) calloc (N,sizeof(float));


	float *in_Dipole ;
	in_Dipole =(float*) calloc (N,sizeof(float)); 
	float *in_Dipole_Variation;
	in_Dipole_Variation=(float*) calloc (N,sizeof(float)); 

	float *out_Poynting;
	out_Poynting=(float*) calloc (N,sizeof(float)); 

	float *out_power;			
	out_power =(float*) calloc (N,sizeof(float));

	
	//////////////////Reading input file ///////////////////////////////
	ifstream input;
	input.open(argv[2]);

string input_name;
input_name=argv[2];

	cout << "Reading input file";
	
	if(input.is_open())
	{
		int i=0;
		while(input.good())
		{	
			// Reading input file 
			input >> Time[i] >>in_net_Q[i] >>in_avg_I[i]>>in_sto_E[i] >> in_P_Ohm[i] >> in_Dipole[i] >> in_Dipole_Variation[i];
			i++;	
		}
	}
	cout << "\tDONE!"<<"\tS="<<S<<endl;
	float dt = Time[2]- Time[1];
	cout << "dt " << dt <<endl;
	cout << "time max " << Time[N-1] <<endl;	
	////////////////////////////////////////////////////////////////////

	ofstream logfile;
	logfile.open("electronicAnalysis.log",std::ios_base::app);
	time_t time_raw;
	struct tm * time_info;
	time (&time_raw);
	time_info = localtime (&time_raw);
	char time_stamp [80];
	strftime (time_stamp,80,"%F %H:%M:%S\n",time_info);
	logfile << "\n#Simulation @ " << time_stamp <<endl;
	logfile << "\n#S value \n" << S <<"\n";
	
	string electrofile = "analysis_" + input_name;
	ofstream data_elec;
	data_elec.open(electrofile);
	data_elec << scientific; 
	
	
	ConvolveGauss(1, 50, 65.0, in_sto_E, out_power, N);
	ConvolveGauss(1, 50, 65.0, in_Dipole_Variation, out_Poynting, N);
	//Derivative1D(N, dt,in_sto_E,  out_power2 );

	
	float avg_Power;
	avg_Power = SignalAverage(N, dt, out_power);	
	logfile << "Average power " << avg_Power <<endl;
	
	for(int i=0;i<N;i++){
		data_elec << Time[i]<<"\t"<<out_power[i]/dt<<"\t"<< out_Poynting[i]/dt <<"\n";
	}	
	cout << "\033[1;7;5;33m Program Running \033[0m"<<endl;
	cout << "\033[1A\033[2K\033[1;32mDONE!\033[0m\n";
	cout<<"═══════════════════════════════════════════════════════════════════════════" <<endl;

	free(Time);
	free(in_net_Q);
	free(in_sto_E);
	free(out_power);
	free(in_Dipole);
	free(in_Dipole_Variation);
	free(out_Poynting);

	logfile.close();
	data_elec.close();
	input.close();
	return 0;
}
