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


#include "Tethys1DLib.h"

using namespace std;




int main(int argc, char **argv){
	cout<<"\n" ;
	cout<<"╔═════════════════════════════════════════════════════════════════════════╗\n";
	cout<<"║ \033[1m Time Series analysis for TETHYS                                       \033[0m ║\n";
	cout<<"╚═════════════════════════════════════════════════════════════════════════╝\n";       	

	int N = atoi(argv[1]);	
	
	float S = atof(argv[3]);	
	
	float *Time;			
	Time =(float*) calloc (N,sizeof(float));
	
	float *in_den_0;			
	in_den_0 =(float*) calloc (N,sizeof(float));
	
	float *in_den_L;			
	in_den_L =(float*) calloc (N,sizeof(float));
	
	float *in_vel_0;			
	in_vel_0 =(float*) calloc (N,sizeof(float));
	
	float *in_vel_L;			
	in_vel_L =(float*) calloc (N,sizeof(float));
	
	
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
			input >> Time[i] >> in_den_L[i] >> in_vel_L[i] >> in_den_0[i] >> in_vel_0[i] ;
			i++;	
		}
	}
	cout << "\tDONE!"<<"\tS="<<S<<endl;
	float dt = Time[2]- Time[1];
	cout << "dt " << dt <<endl;
	cout << "time max " << Time[N-1] <<endl;	
	////////////////////////////////////////////////////////////////////

	ofstream logfile;
	logfile.open("Envelope.log",std::ios_base::app);
	time_t time_raw;
	struct tm * time_info;
	time (&time_raw);
	time_info = localtime (&time_raw);
	char time_stamp [80];
	strftime (time_stamp,80,"%F %H:%M:%S\n",time_info);
	logfile << "\n#Simulation @ " << time_stamp ;
	logfile <<"\n"<< input_name<<"\n";
	//logfile << "\n#S value " << S <<"\n";

	float saturation = 0.0;
	float tau = 0.0;
	float error = 0.0;
	
	
	/* Uses the Functions defined in the library to find the extrema  */
	
	/*not needed since density at x=0 is constant by the boundary conditions*/
	//ExtremaFinding(in_den_0, N, S, dt,saturation,tau, error, "extrema_den_0.dat");
	//cout << "Density saturation at x=0: " << saturation << endl;
	//cout << "Time for 99% of saturation: " << tau <<endl;		


//	string str_snd = to_string(S);
//	string nam_post = "S="+str_snd;
	string extrema_den_L_file = "extrema_den_L_" + input_name;
	ExtremaFinding(in_den_L, N, S, dt,saturation,tau, error, extrema_den_L_file);
	logfile << "#Density saturation at x=L:\n" << saturation << endl;
	logfile << "#Time for 99% of saturation:\n" << tau <<endl;		

	//string str_snd = to_string(S);
	//string nam_post = "S="+str_snd;
	string extrema_vel_0_file = "extrema_vel_0_" + input_name;
	ExtremaFinding(in_vel_0, N, S, dt,saturation,tau, error, extrema_vel_0_file);
	logfile << "#Velocity saturation at x=0:\n" << saturation << endl;
	logfile << "#Time for 99% of saturation:\n" << tau <<endl;		


	//string str_snd = to_string(S);
	//string nam_post = "S="+str_snd;
	string extrema_vel_L_file = "extrema_vel_L_" + input_name;
	ExtremaFinding(in_vel_L, N, S, dt,saturation,tau, error, extrema_vel_L_file);	
	logfile << "#Velocity saturation at x=L:\n" << saturation << endl;
	logfile << "#Time for 99% of saturation:\n" << tau <<endl;		
		
	/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	
	cout << "\033[1;7;5;33m Program Running \033[0m"<<endl;
	cout << "\033[1A\033[2K\033[1;32mDONE!\033[0m\n";
	cout<<"═══════════════════════════════════════════════════════════════════════════" <<endl;

	free(Time);
	free(in_vel_0);
	free(in_vel_L);
	free(in_den_0);
	free(in_den_L);
	logfile.close();
	input.close();
	return 0;
}
