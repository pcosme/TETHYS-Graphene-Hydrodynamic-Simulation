#ifndef TETHYS2DLIB_H
#define TETHYS2DLIB_H

#include <H5Cpp.h>
#include "TethysLib.h"
#include "TethysMathLib.h"

using namespace H5;


class Fluid2D : public TethysBase
{
	protected:
		float * vel_snd_arr;    // array for saving the (potentially varying) S(x,y) function at main grid
		float * vel_snd_arr_mid;// array for saving the (potentially varying) S(x,y) function at auxiliary grid
		float * den_mid ;       // mid or auxiliary grids defined with (Nx-1)*(Ny-1) size
		float * flxX_mid ;
		float * flxY_mid ;

		float * lap_flxX ;      // mass density flux laplacian component x
		float * lap_flxY ;      // mass density flux laplacian component y
		std::ofstream data_preview; // file stream for simplified .dat file output
		int snapshot_per_period = 20;
		int snapshot_step = 1;
	public :
	float * den_old ;       // auxiliary grids for t=n-1 defined with (Nx)*(Ny) size
	float * flxX_old ;
	float * flxY_old ;
	float * den_new ;       // auxiliary grids for t=n-1 defined with (Nx)*(Ny) size
	float * flxX_new ;
	float * flxY_new ;

		float * Den ;       // number density
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
		virtual void SetSimulationTime();   // Finds and set the appropriate simulation time that is 1) Longer than the saturation time 2) Contains enough oscillation periods in the saturated region
		void InitialCondRand();             // Initial condition, zero velocity and constant density with 0.5% white noise
		void InitialCondTest();             // Initial condition for testing and debugging
		virtual void CflCondition();    // Calculates dx and imposes Courant–Friedrichs–Lewy condition to dt
		void Richtmyer();                   // Central Algorithm for solving the hyperbolic conservation law
		virtual float DensityFluxX(__attribute__((unused)) float n, float flx_x, __attribute__((unused)) float flx_y, __attribute__((unused)) float mass, __attribute__((unused)) float s); // density equation (continuity equation) conserved flux X component
		virtual float DensityFluxY(__attribute__((unused)) float n, __attribute__((unused)) float flx_x, float flx_y, __attribute__((unused)) float mass, __attribute__((unused)) float s); // density equation (continuity equation) conserved flux Y component
		virtual float DensitySource(__attribute__((unused)) float n,__attribute__((unused)) float flx_x,__attribute__((unused)) float flx_y,__attribute__((unused)) float mass,__attribute__((unused)) float s);
		virtual float MassFluxXFluxX(float n, float flx_x, float flx_y,__attribute__((unused)) float mass, float s); // velocity X component equation (momentum equation) conserved flux X component
		virtual float MassFluxXFluxY(float n, float flx_x, float flx_y,__attribute__((unused)) float mass, float s); // velocity X component equation (momentum equation) conserved flux Y component
		virtual float MassFluxXSource(__attribute__((unused))float n, float flx_x,__attribute__((unused)) float flx_y,__attribute__((unused)) float mass,__attribute__((unused)) float s);
		virtual float MassFluxYFluxX(float n, float flx_x, float flx_y,__attribute__((unused)) float mass, float s); // velocity Y component equation (momentum equation) conserved flux X component
		virtual float MassFluxYFluxY(float n, float flx_x, float flx_y,__attribute__((unused)) float mass, float s); // velocity Y component equation (momentum equation) conserved flux Y component
		virtual float MassFluxYSource(__attribute__((unused))float n, float flx_x,__attribute__((unused)) float flx_y,__attribute__((unused)) float mass,__attribute__((unused)) float s);
		virtual void MassFluxToVelocity(); // Converts the mass flux density p=mnv to velocity
		void VelocityToCurrent(); // Converts the mass flux density p=mnv to velocity
		void VelocityLaplacianFtcs();
		void VelocityLaplacianDuFortFrankel();

		void TimeUpdate();

		void CreateFluidFile();     // create and open the simplified .dat file output
		void WriteFluidFile(float t) ; // writes the line of time t on the simplified .dat file output
		void SaveSnapShot();
		void ReadSnapShot(const H5std_string &snap_name);
		void SaveSound();
		int GetSnapshotStep() const;
		void SetSnapshotStep() ;
		void SetSnapshotStep(int step) ;
		int GetSnapshotFreq() const;
};

class GrapheneFluid2D : public Fluid2D{
	public :
		explicit GrapheneFluid2D(SetUpParameters &input_parameters);
		~GrapheneFluid2D();
		void CflCondition() override;
		void SetSimulationTime() override;
		void MassFluxToVelocity() override; // Converts the mass density flux back to velocity, in graphene  v = p n^{-3/2}
		/*Override fluxes and sources to specifics of graphene physics*/
		float DensityFluxX(float n, float flx_x, float flx_y,float mass, float s) override;
		float DensityFluxY(float n, float flx_x, float flx_y,float mass, float s) override;
		float DensitySource(float n, float flx_x, float flx_y, float mass, float s)override;
		float MassFluxXFluxX(float n, float flx_x, float flx_y,float mass, float s) override;
		float MassFluxXFluxY(float n, float flx_x, float flx_y,float mass, float s) override;
		float MassFluxXSource(float n, float flx_x, float flx_y, float mass, float s)override;
		float MassFluxYFluxX(float n, float flx_x, float flx_y,float mass, float s) override;
		float MassFluxYFluxY(float n, float flx_x, float flx_y,float mass, float s) override;
		float MassFluxYSource(float n, float flx_x, float flx_y, float mass, float s)override;

		void MagneticSourceSemiAnalytic(); // Semi analytic method for the magnetic interaction
		//void MagneticSourceFtcs();  // Forward Time Centered Space method for the magnetic interaction
		void ParabolicOperatorFtcs();       // Forward Time Centered Space method for the viscous terms
		void ParabolicOperatorDuFortFrankel();       // DuFort Frankel method for the viscous terms
};


#endif
