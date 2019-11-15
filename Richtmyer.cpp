#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <string>
#include <iomanip>   

#include "dyakonovshur.h"


using namespace std;


//<<<<<<< add_source_term
float DensityFlux(float den,float vel,float vel_snd,float vel_fer);

float VelocityFlux(float den,float vel,float vel_snd,float vel_fer);

float DensitySource(float den,float vel,float vel_snd,float vel_fer);

float VelocitySource(float den,float vel,float vel_snd,float vel_fer,float col_freq);



int main(int argc, char **argv){

    BannerDisplay();

	/*......TIME stamp for the logfile................................*/
	ofstream logfile;
	logfile.open("Simulation.log",std::ios_base::app);
	time_t time_raw;
	struct tm * time_info;
	time (&time_raw);
	time_info = localtime (&time_raw);
	logfile << "\n#Simulation @ " << asctime(time_info) ;
	/*................................................................*/
	
	int Nx=201; 							// number of spatial points
	float t=0.0,leng=1.0;					// time variable and spatial Length
	float dx;								// spatial discretisation
	float dt;								// time step
	float vel_snd;						    // Sound speed
	float vel_fer;							// Fermi velocity
	float col_freq; 								// mean free path in units of GFET length
	float *den;							 	//density field
	den =(float*) calloc (Nx,sizeof(float));
	float *den_mid;							//density auxiliary vector for midpoint calculation 
	den_mid = (float*) calloc (Nx-1,sizeof(float));
	float *vel;								//velocity field
 	vel = (float*) calloc (Nx,sizeof(float));
	float *vel_mid;							//velocity auxiliary vector for midpoint calculation 
	vel_mid = (float*) calloc (Nx-1,sizeof(float));
	
	float *den_cor;							//density corrected after average filter 
	den_cor = (float*) calloc (Nx,sizeof(float));
	float *vel_cor;							//velocity corrected after average filter 
	vel_cor = (float*) calloc (Nx,sizeof(float));
 	float *cur_cor;							//current density (n*v) corrected after average filter 
	cur_cor = (float*) calloc (Nx,sizeof(float));
 	
 	
 	
 	int data_save_mode=0;
	
	if(argc==5){
		vel_snd = atof(argv[1]);
		vel_fer = atof(argv[2]);
		col_freq     = atof(argv[3]); 
		data_save_mode = atoi(argv[4]);	// full data or light save option
		}
	else{
		cout << "Define S value: ";
		cin >> vel_snd;
		cout << "Define Vf value: ";
		cin >> vel_fer;
		cout << "Define mean free path value: ";
		cin >> col_freq;
		cout << "Define data_save_mode value (0-> light save | 1-> full data): ";
		cin >> data_save_mode;
		}
	
	
	// NEEDS REVIEWING
	/*......CFL routine to determine dt...............................*/	
	dx = leng / ( float ) ( Nx - 1 );
	if(vel_snd<5){
		dt = dx / (5*vel_snd);
	}
	else{
		if(vel_snd>8 && vel_snd<10){
			dt = dx / (30+3*vel_snd);
		}
		else{
			dt = dx / (20+2*vel_snd);		
			//dt = dx / (5+1.5*vel_snd);		
		}
	}
	/*................................................................*/
	
	
	/*.........Fixed or variable vel_snd value..............................*/
	float *arr_snd;							
	arr_snd =(float*) calloc (Nx,sizeof(float));	
	for(int i = 0; i<Nx  ;i++){
		//arr_snd[i]= vel_snd - 0.15*vel_snd*( dx*i- floor(dx*i) );
		arr_snd[i]=vel_snd;
	}
	/*................................................................*/


	/*.........Output files and streams...............................*/
	
	string str_snd = to_string(vel_snd);
	str_snd.erase(str_snd.end()-4,str_snd.end());
	string str_fer = to_string(vel_fer);
	str_fer.erase(str_fer.end()-4,str_fer.end());
	string str_col_freq = to_string(col_freq);
	str_col_freq.erase(str_col_freq.end()-4,str_col_freq.end());
	
	string nam_post = "S="+str_snd+"vF="+str_fer+"l="+str_col_freq;
		
	// density(x,t)
	string densityfile = "density_" + nam_post + ".dat" ;
	ofstream data_density;
	data_density.open (densityfile);
	data_density << fixed ;
	data_density << setprecision(6);
	// velocity(x,t)	
	string velocityfile = "velocity_" + nam_post + ".dat" ;
	ofstream data_velocity;
	data_velocity.open (velocityfile);
	data_velocity << fixed ;
	data_velocity << setprecision(6);
	// current(x,t)	
	string currentfile = "current_" + nam_post + ".dat" ;
	ofstream data_current;
	data_current.open (currentfile);
	data_current << fixed ;
	data_current << setprecision(6);	

	// time density(L,t)-1=U(L,t) current(0,t) electric_dipole_moment(t)  derivative_electric_dipole_moment(t)
	string electrofile = "electro_" + nam_post + ".dat" ;
	ofstream data_electro;
	data_electro.open (electrofile);
	data_electro << scientific; 
	// time density(L,t) velocity(L,t) density(0,t) velocity(0,t)
	string slicefile = "slice_" + nam_post + ".dat" ;
	ofstream data_slice;
	data_slice.open (slicefile);
	data_slice << scientific; 
	/*................................................................*/

	
	float T_max=10.0;
	
	WellcomeScreen(vel_snd, vel_fer, col_freq,dt,dx,T_max);

/*	cout << "Sound speed S/v0\t"<< vel_snd <<endl;
	cout << "Fermi velocity vF/v0\t"<< vel_fer <<endl;
	if ( PhaseVel(vel_snd, vel_fer) < vel_fer){
		cout << "Phase velocity\t" << PhaseVel(vel_snd, vel_fer)<<"\t WARNING plasmon wave in critical damping region"<<endl;
	}else{
		cout << "Phase velocity\t" << PhaseVel(vel_snd, vel_fer)<<endl;
	}
	cout << "Collision frequency \t"<< col_freq <<endl;
	cout <<"dt= "<<dt<<"\tdx= "<<dx<<endl;
	cout << "Predicted w'= "<< RealFreq(vel_snd,vel_fer,col_freq,1) << "\t1/w'= "<< 1.0/RealFreq(vel_snd,vel_fer,col_freq,1)  << endl;
	cout << "Predicted w''= "<< ImagFreq(vel_snd,vel_fer,col_freq) <<"\t1/w''= "<< 1.0/ImagFreq(vel_snd,vel_fer,col_freq) <<endl;
	
	logfile << "#vel_snd \t vel_fer \t col_freq \t dt \t dx \t w' \t w'' " << endl;
	logfile << vel_snd <<"\t"<<vel_fer<< "\t"<< col_freq<<"\t"<< dt <<"\t"<< dx <<"\t"<< RealFreq(vel_snd,vel_fer,col_freq,1) <<"\t"<< ImagFreq(vel_snd,vel_fer,col_freq) ;
*/	

	
	
	
	
	////////////////////////////////////////////////////////////////////
	// Initialization	
	InitialCondRand(Nx, dx, den, vel);
	BoundaryCond(3, Nx, den, vel);
	////////////////////////////////////////////////////////////////////
	
	if(data_save_mode){
		for(int i = 0; i<Nx  ;i++)
		{
			data_density   <<  den[i] <<"\t";
			data_current   <<  vel[i]*den[i] <<"\t";
			data_velocity  <<  vel[i] <<"\t";
		}
	}
	
	cout << "Running"<<endl;
	
	int time_step=0;
	
	while(t<=T_max && isfinite(vel[(Nx-1)/2]))
	{	
		++time_step;
		t += dt;
		
		if(data_save_mode && time_step % 35 == 0 ){
			data_density  << "\n";
			data_current  << "\n";
			data_velocity << "\n";
		}
		
		//
		//  Half step calculate density and velocity at time k+0.5 at the spatial midpoints
		//
		for ( int i = 0; i < Nx - 1; i++ )
		{
			den_mid[i] = 0.5*( den[i] + den[i+1] )
				- ( 0.5*dt/dx ) * ( DensityFlux(den[i+1],vel[i+1],arr_snd[i], vel_fer) - DensityFlux(den[i],vel[i],arr_snd[i], vel_fer) ) 
				+ ( 0.5*dt    ) * DensitySource(0.5*(den[i]+den[i+1]),0.5*(vel[i]+vel[i+1]),arr_snd[i], vel_fer) ;
			vel_mid[i] = 0.5*( vel[i] + vel[i+1] )
				- ( 0.5*dt/dx ) * ( VelocityFlux(den[i+1],vel[i+1],arr_snd[i], vel_fer) - VelocityFlux(den[i],vel[i],arr_snd[i], vel_fer) ) 
				+ ( 0.5*dt    ) * VelocitySource(0.5*(den[i]+den[i+1]),0.5*(vel[i]+vel[i+1]),arr_snd[i], vel_fer, col_freq) ;
		}
		//
		// Remaining step 
		//
		for ( int i = 1; i < Nx - 1; i++ )
		{
			den[i] = den[i] - (dt/dx) * ( DensityFlux(den_mid[i],vel_mid[i],arr_snd[i], vel_fer) - DensityFlux(den_mid[i-1],vel_mid[i-1],arr_snd[i], vel_fer) )
							+  dt * DensitySource(den[i],vel[i],arr_snd[i], vel_fer);
			vel[i] = vel[i] - (dt/dx) * ( VelocityFlux(den_mid[i],vel_mid[i],arr_snd[i], vel_fer) - VelocityFlux(den_mid[i-1],vel_mid[i-1],arr_snd[i], vel_fer) )
							+  dt * VelocitySource(den[i],vel[i],arr_snd[i], vel_fer, col_freq);
		}
		
		// Impose boundary conditions
		BoundaryCond(3, Nx, den, vel);
		
		// Applying average filters for smoothing 
		AverageFilter( den ,den_cor, Nx , 2);	
		AverageFilter( vel ,vel_cor, Nx , 2);
		
		// calculate current density already smoothed			
		for ( int i = 0; i < Nx; i++ )
		{	
			cur_cor[i] = 	vel_cor[i]*den_cor[i];	
			if(data_save_mode){
			//Record full data
				if(time_step % 35 == 0){
					data_density   << den_cor[i] <<"\t";
					data_current   << cur_cor[i] <<"\t";
					data_velocity  << vel_cor[i] <<"\t";
				}
			}
		}
	
		//Record end points
		data_slice <<t<<"\t"<< den_cor[Nx-1] <<"\t"<< vel_cor[Nx-1] <<"\t"<< den_cor[0] <<"\t" << vel_cor[0] <<"\n";
		//Record electric quantities
		data_electro <<t<<"\t"<< den_cor[Nx-1]-1.0 <<"\t"<< den_cor[0]*vel_cor[0] <<"\t"<<  TotalElectricDipole(Nx,dx,den_cor)<<"\t"<<  DtElectricDipole(Nx,dx,cur_cor) <<"\t"<< KineticEnergy(Nx,dx, den_cor, vel_cor)  <<"\n";
	}
	cout << "DONE!" <<endl;
	cout << "*******************************************************"<<endl;

	free(den);
	free(den_mid);
	free(den_cor);
	free(vel);
	free(vel_mid);
	free(vel_cor);
	free(cur_cor);

	data_density.close();
	data_velocity.close();
	data_current.close();
	data_slice.close();
	data_electro.close();
	
	
	
	
	return 0;
}

float DensityFlux(float den,float vel,float vel_snd,float vel_fer){
	float f1;
	
	f1 = den*vel;
	
	return f1;
}


float VelocityFlux(float den,float vel,float vel_snd,float vel_fer){
	float f2;
	f2 = 0.25*vel*vel + vel_fer*vel_fer*0.5*log(den) + 2*vel_snd*vel_snd*sqrt(den); 
	return f2;
}

float DensitySource(float den,float vel,float vel_snd,float vel_fer){
	float Q1=0.0;
return Q1;	
}

float VelocitySource(float den,float vel,float vel_snd,float vel_fer,float col_freq){
	float Q2=0.0;
	
		Q2=-1.0*col_freq*(vel-1);
	
return Q2;
}

