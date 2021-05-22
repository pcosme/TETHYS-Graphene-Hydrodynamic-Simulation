//
// Created by pcosme on 23/12/2020.
//

#ifndef SETUPPARAMETERSLIB_H
#define SETUPPARAMETERSLIB_H
#include <fstream>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <algorithm>
#include <string>
#include <random>
#include <exception>
#include <H5Cpp.h>
#include <omp.h>
/*!
 * @brief Initialization class for the fluid classes.
 *
 * This class allows the initialization of the appropriate physical parameters either from inline arguments, prompt or reading them from an existing HDF5 file
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
		float OddViscosity;
		float CyclotronFrequency;
		float ThermalDiffusivity;
		void ParametersChecking() const; ///< Runs a checking on the physical feasibility of the parameters
		void DefineGeometry(); ///< Set ups the 2D grid dimensions
		void ParametersFromHdf5File(const std::string& hdf5name); ///< Imports the parameters from a saved HDF5 file
        void GetParameters();
};



#endif //SETUPPARAMETERSLIB_H
