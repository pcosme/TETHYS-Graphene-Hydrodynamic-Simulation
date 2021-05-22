/************************************************************************************************\
* Copyright (c) 2020 Pedro Cosme and João Santos                                                 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/

#ifndef FLUID2DLIB_H
#define FLUID2DLIB_H

#include <H5Cpp.h>
#include "TethysBaseLib.h"
#include "TethysMathLib.h"
#include "SetUpParametersLib.h"
#include "Grid2DLib.h"

using namespace H5;


/*!
 * @brief Generalistic fluid class in two dimensions, mainly for testing purposes.
 *
 * The Fluid2D class describes a regular newtonian compressible fluid governed by the usual continuity and Cauchy momemtum equations, where the mass of the fluid element is constant.
 * It is maintained mainly for testing, the GrapheneFluid2D class is derived from it and overrides the necessary methods in order to describe the semi-classical electronic fluid.
 * */
class Fluid2D : public TethysBase
{
	protected:
		float * vel_snd_arr;    // array for saving the (potentially varying) S(x,y) function at main grid
		float * vel_snd_arr_mid;// array for saving the (potentially varying) S(x,y) function at auxiliary grid
		float * den_mid ;       // mid or auxiliary grids defined with (Nx-1)*(Ny-1) size
		float * tmp_mid ;
		float * flxX_mid ;
		float * flxY_mid ;

		float * velX_dx;
		float * velX_dy;

		float * velY_dx;
		float * velY_dy;

		float * velX_dx_mid;
		float * velX_dy_mid;

		float * velY_dx_mid;
		float * velY_dy_mid;

		float * lap_flxX ;      // mass density flux laplacian component x
		float * lap_flxY ;      // mass density flux laplacian component y
        float * lap_tmp ;      // temperature laplacian

    std::ofstream data_preview; // file stream for simplified .dat file output
		int snapshot_per_period = 40;
		int snapshot_step = 1;
		void ForwardTimeOperator();


		void VelocityGradient();

		void VelocityGradientMid();

		virtual float DensityToMass(float density);

public :
		float * Den ;       // number density
		float * Tmp ;       // electron temperature
		float * VelX ;      // fluid velocity x component
		float * VelY ;      // fluid velocity y component
		float * FlxX ;      // mass density flux x component
		float * FlxY ;      // mass density flux y component
		float * CurX ;      // current density x component
		float * CurY ;      // current density y component
		explicit Fluid2D(const SetUpParameters &input_parameters);
		~Fluid2D();
		bool Snapshot() const;

		void SetSound();     // Applies the anisotropy to the sound velocity array
		virtual void SetSimulationTime();   ///< Finds and set the appropriate simulation time
		void InitialCondRand();             ///< Initial condition, zero velocity and constant density with 0.5% white noise
		void InitialCondTest();             // Initial condition for testing and debugging
		/*!
		 * @brief Calculates @f$\Delta x@f$ and imposes Courant–Friedrichs–Lewy condition to @f$\Delta t@f$
		 *
		 * The method retrieves the spatial discretization @f$\Delta x = L / N_x @f$ and then sets the time sted as @f$\Delta t = \Delta x /10 @f$.
		 * This in general will suffice for the tests and runs performed solely with the Fluid2D class.
		 */
		virtual void CflCondition();

		/*!
		 * @brief Richtmyer method to solve the hyperbolic fluid equations
		 *
		 * Central algorithm for solving the 2D hyperbolic conservation law equations @cite LeVeque1992.
		 * The behaviour of this solver is controlled by the definitions of the appropriate flux and source functions, cf. the @ref twodrich "descriptive page"
		 *
		 *
		 *
		 * @see DensityFluxX
		 * @see DensityFluxY
		 * @see XMomentumFluxX
		 * @see XMomentumFluxY
		 * @see YMomentumFluxX
		 * @see YMomentumFluxY
		 * @see DensitySource
		 * @see XMomentumSource
		 * @see YMomentumSource
		 *
		 */
		void Richtmyer();                   // Central Algorithm for solving the hyperbolic conservation law


//		virtual float DensityFluxX( float n, float flx_x,  float flx_y,  float mass,  float s); ///< density equation (continuity equation) conserved flux X component
//		virtual float DensityFluxY( float n,  float flx_x, float flx_y,  float mass,  float s); ///< density equation (continuity equation) conserved1 flux Y component
//		virtual float XMomentumFluxX(float n, float flx_x, float flx_y, float mass, float s); ///< velocity X component equation (momentum equation) conserved flux X component
//		virtual float XMomentumFluxY(float n, float flx_x, float flx_y, float mass, float s); ///< velocity X component equation (momentum equation) conserved flux Y component
//		virtual float YMomentumFluxX(float n, float flx_x, float flx_y, float mass, float s); ///< velocity Y component equation (momentum equation) conserved flux X component
//		virtual float YMomentumFluxY(float n, float flx_x, float flx_y, float mass, float s); ///< velocity Y component equation (momentum equation) conserved flux Y component

	virtual float DensitySource( float n, float flx_x, float flx_y, float mass, float s); ///< density equation (continuity equation) source term
	virtual float XMomentumSource(float n, float flx_x, float flx_y, float mass, float s); ///< velocity X component equation (momentum equation) source term
	virtual float YMomentumSource(float n, float flx_x, float flx_y, float mass, float s); ///< velocity y component equation (momentum equation) source term


    virtual float DensityFluxX(GridPoint p, char side ); ///< density equation (continuity equation) conserved flux X component
    virtual float DensityFluxY(GridPoint p, char side ); ///< density equation (continuity equation) conserved1 flux Y component

    virtual float XMomentumFluxX(GridPoint p, char side ); ///< velocity X component equation (momentum equation) conserved flux X component
    virtual float XMomentumFluxY(GridPoint p, char side ); ///< velocity X component equation (momentum equation) conserved flux Y component

    virtual float YMomentumFluxX(GridPoint p, char side ); ///< velocity Y component equation (momentum equation) conserved flux X component
    virtual float YMomentumFluxY(GridPoint p, char side ); ///< velocity Y component equation (momentum equation) conserved flux Y component

    float TemperatureFluxX(GridPoint p, char side );
    float TemperatureFluxY(GridPoint p, char side );



		/*!
		* @brief Converts the mass density flux to velocity on the entire simulation grid.
		*
		* Since the mass of the fluid element is constant one needs only to perform the transformation
		@f[ \vec{v} = \frac{\vec{p}}{n} @f]
		* */
		virtual void MassFluxToVelocity(); // Converts the mass flux density p=mnv to velocity

		/*!
		* @brief Converts velocity field to current density on the entire simulation grid.
		*
		* The method simply performs
		@f[ \vec{j} = \vec{v}n @f]
		* */
		void VelocityToCurrent(); // Converts the mass flux density p=mnv to velocity


		void CreateFluidFile();     ///< creates and opens the simplified .dat file output

		/*!
		 * @brief Writes the line of time t on the simplified .dat file output
		 *
		 * As a way to easily access the main results of the simulation the simplified .dat file stores the following quantities by columns:
		 * -# Time @f$t@f$
		 * -# Density at drain contact @f$n(x=L)@f$
		 * -# Mass flux along x at drain contact  @f$p_x(x=L)@f$
		 * -# Density at source contact @f$n_x(x=0) @f$
		 * -# Mass flux along x at source contact @f$p_x(x=0) @f$
		 * */
		void WriteFluidFile(float t) ; // writes the line of time t on the simplified .dat file output

		/*!
		 * @brief Saves the current snapshot on HDF5 file
		 *
		 *
		 * */
		void SaveSnapShot();

		/*!
		 * @brief Imports snapshot to a Fluid2D class
		 *
		 * @param snap_name Name of the snapshot to import
		 * */
		void ReadSnapShot(const H5std_string &snap_name);

		/*!
		 * @brief Saves the grid of the spatial distribution of the parameter S at the output HDF5 file
		 *
		 * */
		void SaveSound();
		int GetSnapshotStep() const; ///< Returns the number of the present snapshot @return snapshot_step order number of the snapshot
		int GetSnapshotFreq() const; ///< Returns the number of snapshots per period to record  @return snapshot_per_period number of the snapshots per period



		/*!
		 * @brief Calculates the velocity Laplacians for the FTCS method
		 *
		 * Them method computes both laplacians @f$ \nabla^2 v_x @f$ and @f$ \nabla^2 v_y @f$ using a five points second order stencil
		 * @f[\nabla^2 u \approx \frac{1}{\Delta x^2} ( -4u_{i,j}+u_{i-1,j}+u_{i+1,j}+u_{i,j-1}+u_{i,j+1} )  @f]
		 * @see ParabolicOperatorFtcs
 		 * */
		void VelocityLaplacianFtcs();
		/*!
		* @brief Calculates the velocity Laplacians for the Weighted (1,9) method
		*
		* The method computes both laplacians @f$ \nabla^2 v_x @f$ and @f$ \nabla^2 v_y @f$ using a nine points fourth order stencil. More details can be found, instance in:
		* <a href="https://digital.library.adelaide.edu.au/dspace/bitstream/2440/18689/2/02whole.pdf">this work</a>
		*
		* @see ParabolicOperatorWeightedExplicit19()
		* */
		float Laplacian19(GridPoint p, float *input_ptr, float constant);

		void VelocityLaplacianWeighted19();
		/*!
		* @brief Forward Time Centered Space method for the viscous terms
		*
		* @see VelocityLaplacianFtcs()
		* */
		void ParabolicOperatorFtcs();       // Forward Time Centered Space method for the viscous terms
		/*!
		* @brief Forward Time Weighted (1,9) space method
		*
		* @see VelocityLaplacianWeighted19()
		* */
		void ParabolicOperatorWeightedExplicit19(); // Forward Time Centered Space method for the viscous terms
        void TemperatureLaplacianWeighted19();
};


#endif
