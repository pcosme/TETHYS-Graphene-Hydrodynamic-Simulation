#ifndef TETHYSLIB_H
#define TETHYSLIB_H

#include <fstream>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <string>
#include <random>
#include <exception>

#include <H5Cpp.h>
using namespace std;
using namespace H5;

#ifndef MAT_PI
#	define MAT_PI 3.14159265358979f
#endif

#ifndef MAT_EULER
#	define MAT_EULER 2.71828182845905f
#endif

const FloatType      HDF5FLOAT(PredType::NATIVE_FLOAT);
const IntType        HDF5INT(PredType::NATIVE_INT);


/*
 * Struct to pass the initialization
 * */
class SetUpParameters {
	public:
		SetUpParameters();
		SetUpParameters(int argc, char ** argv);
		SetUpParameters(float sound, float fermi, float coll, float visco, float cyclo, int mode, float aspect);
		~SetUpParameters() = default;
		int SaveMode;
		int SizeX;
		int SizeY;
		float Length=1.0f;
		float Width=1.0f;
		float AspectRatio=1.0f;
		float SoundVelocity;
		float FermiVelocity;
		float CollisionFrequency;
		float ShearViscosity;
		float CyclotronFrequency;
		void ExceptionsChecking() const;
		void DefineGeometry();
		void ParametersFromHdf5File(const std::string& hdf5name);
};


/* Base class from which the fluid and graphene fluid classes are derived
 * The setting of the dimensions and creation of the HDF5 file and its attributes is controlled from this class
 * */
class TethysBase {
	protected:
		int   Nx ;          // Simulation region (dataset) dimensions
		int   Ny ;
		int   RANK;         // Either 1 or 2 indicating the dimension of the problem
		float dx=1.0f;      // space discretization along x
		float dy=1.0f;      // space discretization along y
		float dt=1.0f;      // temporal discretization. will later be redifined by the CFL condition
		float lengX=1.0f;   // physical length along x. dx will be redifined as lengX/Nx
		float lengY=1.0f;   // physical length along y. dy will be redifined as lengY/Ny
		float vel_snd =50.0f;   // sound velocity parameter
		float vel_fer =10.0f;
		float cyc_freq =0.0f;
		float kin_vis =0.0f;    // kinetic shear viscosity parameter
		float col_freq =0.0f;   // colision frequency parameter
		std::string file_infix; // base name for the output files
		float Tmax=10;          // total time of simulation

	public:
		TethysBase(int size_nx, int size_ny, int dimensions); // class constructor initializes Nx, Ny, RANK and file_infix
		~TethysBase();  // class destructor if the the flag HDF5fileOpen=TRUE it deletes the dataspaces and hdf5 files

		static bool HDF5fileOpen;
		static int TimeStepCounter;
		static float TimeStamp;

		H5File* Hdf5File ;  // hdf5 file handler
		Group* GrpDat ;     // group for the simulated data: Attributes; Density; Velocity X; Velocity Y
		Group* GrpDen ;     // group for ALL Density snapshots
		Group* GrpVelX ;    // group for ALL Velocity X snapshots
		Group* GrpVelY ;    // group for ALL Velocity X snapshots
		DataSpace* DataspaceVelSnd; // dataspace for the sound anisotropy
		DataSpace* DataspaceVelSndMid; // dataspace for the sound anisotropy
		DataSpace* DataspaceDen;    // dataspace for EACH Density snapshots
		DataSpace* DataspaceVelX;   // dataspace for EACH Velocity X snapshots
		DataSpace* DataspaceVelY;   // dataspace for EACH Velocity Y snapshots

		void SetTmax(float x);      // setter method for the total simulation time cf. GrapheneFluid2D::SetSimulationTime()
		float GetTmax() const;      // getter method for the total simulation time
		int SizeX() const;          // getter method for the number of simulation points along x
		int SizeY() const;          // getter method for the number of simulation points along y
		int Rank() const;           // getter method for the system dimensionality

		void SetVelSnd(float x);    // setter method for nominal S value
		void SetKinVis(float x);    // setter method for kinetic shear viscosity
		void SetColFreq(float x);   // setter method for collision frequency
		void SetVelFer(float x);        // setter method for Fermi Velocity
		void SetCycFreq(float x);        // setter method for cyclotron frequency
		void SetDx(float x);        // setter method for spatial step x
		void SetDy(float x);        // setter method for spatial step y
		void SetDt(float x);        // setter method for temporal step
		void SetLengthX(float x);   // setter method for total length along x
		void SetLengthY(float x);   // setter method for total length along y

		float GetVelSnd() const;    // getter method for nominal S value
		float GetKinVis() const;    // getter method for kinetic shear viscosity
		float GetColFreq() const;   // getter method for collision frequency
		float GetVelFer() const;        // getter method for Fermi Velocity
		float GetCycFreq() const;       // getter method for cyclotron frequency
		float GetDx() const;        // getter method for spatial discretization x
		float GetDy() const;        // getter method for spatial discretization y
		float GetDt() const;        // getter method for time discretization
		float GetLengthX() const;   // getter method for total length along x
		float GetLengthY() const;   // getter method for total length along y

		float ImagFreq() const;
		float PhaseVel() const;
		float RealFreq() const;


		std::string GetInfix() const;   // getter method for file name infix


		void CreateHdf5File();          // creates the HDF5 files with the necessary structure
		void OpenHdf5File(const std::string& hdf5name); // opens an existing HDF5 file with the necessary structure
		void CloseHdf5File() const;           // closes the HDF5 file
		void WriteAttributes();          // saves the simulation attributes (either physical and simulation parameters)

		void BannerDisplay() const; // launches the initial ASCII art banner
		void WelcomeScreen() const; //launches screen with the relevant info
};
#endif

