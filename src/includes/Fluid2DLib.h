/************************************************************************************************\
* 2020 Pedro Cosme , João Santos, Ivan Figueiredom, João Rebelo, Diogo Simões                    *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/


/*!@file
 * @brief Header file for 1D general fluid
 */


#ifndef FLUID2DLIB_H
#define FLUID2DLIB_H

#include <functional>
#include <H5Cpp.h>
#include "includes/TethysBaseLib.h"
#include "includes/TethysMathLib.h"
#include "includes/SetUpParametersLib.h"
#include "includes/GridLib.h"
#include "includes/StateVec2DLib.h"

using namespace H5;

/*!
 * @brief Generalistic fluid class in two dimensions, mainly for testing purposes.
 *
 * The Fluid2D class describes a regular newtonian compressible fluid governed by the usual continuity and Cauchy momemtum equations, where the mass of the fluid element is constant.
 * It is maintained mainly for testing, the GrapheneFluid2D class is derived from it and overrides the necessary methods in order to describe the semi-classical electronic fluid.
 * */
	class Fluid2D : public TethysBase {
	protected:

		std::ofstream data_preview; // file stream for simplified .dat file output
		int snapshot_per_period = 40;
		int snapshot_step = 1;

		virtual void ForwardTimeOperator(); ///< Time evolution for the FTCS method employed for the parabolic operators.
		virtual void ForwardTimeOperator(char field); ///< Time evolution for the FTCS method employed for the parabolic operators.


		virtual float DensityToMass(float density);

		virtual void ChooseGridPointers(const string &grid);

		float SideAverage(const float *input_array, GridPoint2D p, char side);
		StateVec2D SideAverage(const StateVec2D *input_array, GridPoint2D p, char side);
		StateVec2D *ptr_StateVec = nullptr;

		virtual void RichtmyerStep1();
		virtual void RichtmyerStep2();

	public :

		virtual void CopyFields();

		StateVec2D * Umain;
		StateVec2D * Umid;

		float *vel_snd_arr;    // array for saving the (potentially varying) S(x,y) function at main grid

		float *Den;       // number density
		float *Tmp;       // electron temperature
		float *VelX;      // fluid velocity x component
		float *VelY;      // fluid velocity y component

		float *CurX;      // current density x component
		float *CurY;      // current density y component

		explicit Fluid2D(const SetUpParameters &input_parameters);
		~Fluid2D();

		bool Snapshot() const;



		void SetSound();     ///< sets a constant sound velocity
		/*!
		* @brief 2D Function of local anisotropy of S
		*
		* Function to implement the spatial variation of the sound velocity S(x,y) in 2D
		* corresponding to a variation of substrate permittivity or even the description of a multi gated system.
		*
		* @param std::function<float(float,float)> real function of real arguments f(x,y)
		*
		* */
		void SetSound(std::function<float(float,float)>);     ///< Applies the anisotropy (in the cases there is one) to the sound velocity array
		virtual void SetSimulationTime();   ///< Finds and set the appropriate simulation time

		virtual void InitialCondRand();             ///< Initial condition, zero velocity and constant density with 0.5% white noise
		void InitialCondWave();
		void InitialCondTest();             ///< Initial condition for testing and debugging
		void InitialCondGeneral(function<float(float, float)> fden, function<float(float, float)> fvx,
		                        function<float(float, float)> fvy);

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
		 * @see EleDensityFluxX
		 * @see EleDensityFluxY
		 * @see EleXMomentumFluxX
		 * @see EleXMomentumFluxY
		 * @see EleYMomentumFluxX
		 * @see EleYMomentumFluxY
		 * @see EleDensitySource
		 * @see EleXMomentumSource
		 * @see EleYMomentumSource
		 *
		 */
	    void Richtmyer();                   // Central Algorithm for solving the hyperbolic conservation law

		virtual float DensitySource(StateVec2D U); ///< density equation (continuity equation) source term
		virtual float XMomentumSource(StateVec2D U); ///< velocity X component equation (momentum equation) source term
		virtual float YMomentumSource(StateVec2D U); ///< velocity y component equation (momentum equation) source term
		virtual float TemperatureSource(StateVec2D U); ///< density equation (continuity equation) source term

		virtual float DensityFluxX(StateVec2D U);///< density equation (continuity equation) conserved flux X component
		virtual float DensityFluxY(StateVec2D U); ///< density equation (continuity equation) conserved1 flux Y component
		virtual float XMomentumFluxX(StateVec2D U); ///< velocity X component equation (momentum equation) conserved flux X component
		virtual float XMomentumFluxY(StateVec2D U); ///< velocity X component equation (momentum equation) conserved flux Y component
		virtual float YMomentumFluxX(StateVec2D U); ///< velocity Y component equation (momentum equation) conserved flux X component
		virtual float YMomentumFluxY(StateVec2D U); ///< velocity Y component equation (momentum equation) conserved flux Y component
		float TemperatureFluxX(StateVec2D U); ///< Temperature equation (heat equation) conserved flux X component
		float TemperatureFluxY(StateVec2D U); ///< Temperature equation (heat equation) conserved flux Y component


		/*!
		* @brief Converts the mass density flux to velocity on the entire simulation grid.
		*
		* Since the mass of the fluid element is constant one needs only to perform the transformation
		@f[ \vec{v} = \frac{\vec{p}}{n} @f]
		* */
		void MassFluxToVelocity(); // Converts the mass flux density p=mnv to velocity

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
		virtual void WriteFluidFile(float t); // writes the line of time t on the simplified .dat file output

		/*!
		 * @brief Saves the current snapshot on HDF5 file
		 *
		 *
		 * */
		virtual void SaveSnapShot();

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
		void SaveSound(); ///< Records a mesh of eventual sound velocity anisotropy at the HDF5 file
		int
		GetSnapshotStep() const; ///< Returns the number of the present snapshot @return snapshot_step order number of the snapshot
		int
		GetSnapshotFreq() const; ///< Returns the number of snapshots per period to record  @return snapshot_per_period number of the snapshots per period



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
		float Laplacian19(GridPoint2D p, float *input_ptr, float constant);

		virtual void VelocityLaplacianWeighted19();

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
		virtual void ParabolicOperatorWeightedExplicit19(); ///< Forward Time Centered Space method for the diffusive terms
		virtual void ParabolicOperatorWeightedExplicit19(char field); ///< Forward Time Centered Space method for the diffusive terms
		void TemperatureLaplacianWeighted19();



		void VelocityGradient(StateVec2D *Uarray, int size_x, int size_y);

	};


#endif
