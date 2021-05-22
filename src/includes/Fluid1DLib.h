/************************************************************************************************\
* Copyright (c) 2020 Pedro Cosme and João Santos                                                 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/

#ifndef FLUID1DLIB_H
#define FLUID1DLIB_H

#include <H5Cpp.h>
#include "TethysBaseLib.h"
#include "TethysMathLib.h"
#include "SetUpParametersLib.h"

using namespace H5;
using namespace std;


/*!
 * @brief Generalistic fluid class in one dimension, mainly for testing purposes.
 *
 * The Fluid1D class describes a regular newtonian compressible fluid governed by the usual continuity and Cauchy momemtum equations, where the mass of the fluid element is constant.
 * It is maintained mainly for testing, the GrapheneFluid1D class is derived from it and overrides the necessary methods in order to describe the semi-classical electronic fluid.
 * */
class Fluid1D : public TethysBase{
	protected:
		float * vel_snd_arr;        // array for saving the (potentially varying) S(x) function
		float * den_mid ;           // number density at midpoint
		float * vel_mid ;           // velocity at midpoint
		float * grad_vel_mid;       // velocity gradient at mid point for the viscous case
		std::ofstream data_preview; // file stream for simplified .dat file output
		int snapshot_per_period = 10;
		int snapshot_step = 1;

	public :
		float * Den ;       // number density
		float * Vel ;       // fluid velocity
		float * GradVel;    // fluid velocity gradient
		float * Cur ;       // current density (density times velocity)
		float * DenCor;     // corrected i.e. smoothed quantities
		float * VelCor ;
		float * CurCor ;
		explicit Fluid1D(const SetUpParameters &input_parameters);
		~Fluid1D();
		bool Snapshot() const;
		void Smooth(int width);     // smoothing moving average filter to obtain the "Cor" version of the quantities
		void SetSimulationTime();   // Finds and set the appropriate simulation time that is 1) Longer than the saturation time 2) Contains enough oscillation periods in the saturated region
		void InitialCondRand();     // Initial condition, zero velocity and constant density with 0.5% white noise
		void InitialCondTest();     // Initial condition for testing and debugging
		void Richtmyer();           // Central Algorithm for solving the hyperbolic conservation law
		void SetSound();            // Applies the anisotropy to the sound velocity array
		virtual void CflCondition();    // Calculates dx and imposes Courant–Friedrichs–Lewy condition to dt
		virtual float DensityFlux(float n,float v, __attribute__((unused)) float s);    // density equation (continuity equation) conserved flux
		virtual float VelocityFlux(float n,float v,float dv, __attribute__((unused)) float s); // velocity equation (momentum equation) conserved flux
		virtual float DensitySource(__attribute__((unused)) float n,  __attribute__((unused)) float v, __attribute__((unused)) float s); // density equation (continuity equation) source term
		virtual float VelocitySource( __attribute__((unused)) float n, __attribute__((unused)) float v, __attribute__((unused)) float s); // velocity equation (momentum equation) source term
		void CreateFluidFile();     // create and open the simplified .dat file output
		void SaveSnapShot(); // saves the all the simulated quantities on the appropriate dataspace of the HDF5 file
		void WriteFluidFile(float t) ; // writes the line of time t on the simplified .dat file output
		int GetSnapshotStep() const;
		int GetSnapshotFreq() const;
};

#endif

