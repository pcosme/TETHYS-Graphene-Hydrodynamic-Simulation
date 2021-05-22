/************************************************************************************************\
* Copyright (c) 2020 Pedro Cosme and João Santos                                                 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/

#ifndef TETHYSBASELIB_H
#define TETHYSBASELIB_H

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
#include <functional>

#include <H5Cpp.h>
#include <omp.h>

using namespace std;
using namespace H5;

#ifndef MAT_PI
#	define MAT_PI 3.14159265358979f
#endif

#ifndef MAT_EULER
#	define MAT_EULER 2.71828182845905f
#endif

#ifndef PHYS_FERMI_CNVC
#	define PHYS_FERMI_CNVC 0.07599088773175f
#endif


const FloatType      HDF5FLOAT(PredType::NATIVE_FLOAT);
const IntType        HDF5INT(PredType::NATIVE_INT);


/*!
 * @brief Base class for the fluid classes
 *
 * This base class, from which the subsequent fluid classes are derived, establish the dimensios of the simulation grid and manage the creation of the HDF5 structures
 *
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
		float odd_vis =0.0f;    // kinetic odd viscosity parameter
		float therm_diff = 0.0f; // thermal diffusivity parameter
		float col_freq =0.0f;   // colision frequency parameter
		std::string file_infix; // base name for the output files
		float Tmax=10;          // total time of simulation

	public:
		/*!
		 Class constructor
		 @param size_nx Number of grid points along x
		 @param size_ny Number of grid points along y
		 @param dimensions Dimensionality of the simulation
		*/
		TethysBase(int size_nx, int size_ny, int dimensions);
		~TethysBase();  // class destructor if the the flag Hdf5FileOpen=TRUE it deletes the dataspaces and hdf5 files

		static bool Hdf5FileOpen;
		static int TimeStepCounter;
		static float TimeStamp;

		H5File* Hdf5File ;  ///< hdf5 file handler
		Group* GrpDat ;     ///< group for the simulated data: Attributes; Density; Velocity X; Velocity Y
		Group* GrpDen ;     ///< group for ALL Density snapshots
		Group* GrpVelX ;    ///< group for ALL Velocity X snapshots
		Group* GrpVelY ;    ///< group for ALL Velocity X snapshots
        Group* GrpTmp ;    ///< group for ALL Temperature snapshots

        DataSpace* DataspaceVelSnd; ///< dataspace for the sound anisotropy
		DataSpace* DataspaceVelSndMid; ///< dataspace for the sound anisotropy
		DataSpace* DataspaceDen;    ///< dataspace for EACH Density snapshots
		DataSpace* DataspaceVelX;   ///< dataspace for EACH Velocity X snapshots
		DataSpace* DataspaceVelY;   ///< dataspace for EACH Velocity Y snapshots
        DataSpace* DataspaceTmp;   ///< dataspace for Temperature

		void SetTmax(float x);      ///< sets  the total simulation time cf. GrapheneFluid2D::SetSimulationTime()
		float GetTmax() const;      ///< @returns   the total simulation time
		int SizeX() const;          ///< @returns   the number of simulation points along x
		int SizeY() const;          ///< @returns   the number of simulation points along y
		int Rank() const;           ///< @returns   the system dimensionality

		void SetVelSnd(float x);    ///< sets  nominal S value
		void SetKinVis(float x);    ///< sets  kinetic shear viscosity
		void SetOddVis(float x);    ///< sets  kinetic odd viscosity
		void SetColFreq(float x);   ///< sets  collision frequency
		void SetThermDiff(float x);   ///< sets  collision frequency
		void SetVelFer(float x);        ///< sets  Fermi Velocity
		void SetCycFreq(float x);        ///< sets  cyclotron frequency
		void SetDx(float x);        ///< sets  spatial step x
		void SetDy(float x);        ///< sets  spatial step y
		void SetDt(float x);        ///< sets  temporal step
		void SetLengthX(float x);   ///< sets  total length along x
		void SetLengthY(float x);   ///< sets  total length along y

		float GetThermDiff() const ;    ///< @returns  collision frequency
		float GetVelSnd() const;    ///< @returns   nominal S value
		float GetKinVis() const;    ///< @returns   kinetic shear viscosity
		float GetOddVis() const;    ///< @returns   kinetic odd viscosity
		float GetColFreq() const;   ///< @returns   collision frequency
		float GetVelFer() const;        ///< @returns   Fermi Velocity
		float GetCycFreq() const;       ///< @returns   cyclotron frequency
		float GetDx() const;        ///< @returns   spatial discretization x
		float GetDy() const;        ///< @returns   spatial discretization y
		float GetDt() const;        ///< @returns   time discretization
		float GetLengthX() const;   ///< @returns   total length along x
		float GetLengthY() const;   ///< @returns   total length along y

		float ImagFreq() const;
		float PhaseVel() const;
		float RealFreq() const;


		std::string GetInfix() const;   ///< @returns   file name infix


		void CreateHdf5File();          ///< creates the HDF5 files with the necessary structure
		void OpenHdf5File(const std::string& hdf5name); ///< opens an existing HDF5 file with the necessary structure @param hdf5name HDF5 file name
		void CloseHdf5File() const;           ///< closes the HDF5 file
		void WriteAttributes();          ///< saves the simulation attributes (either physical and simulation parameters)

		static void BannerDisplay() ; ///< launches the initial ASCII art banner
		void WelcomeScreen() const; ///< launches screen with the relevant info
};
#endif

