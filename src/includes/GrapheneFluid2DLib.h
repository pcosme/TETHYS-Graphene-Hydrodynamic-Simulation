//
// Created by pcosme on 27/12/2020.
//

#ifndef GRAPHENEFLUID2DLIB_H
#define GRAPHENEFLUID2DLIB_H

#include "Fluid2DLib.h"

/*!
 * @brief Graphene electronic fluid class in two dimensions.
 *
 * The GrapheneFluid2D class describes a  fluid governed by the usual continuity and Cauchy momemtum equations, where the mass of the fluid element is *not* constant.
 * It overrides class Fluid2D necessary methods in order to describe the semi-classical electronic fluid.
 * */
class GrapheneFluid2D : public Fluid2D{

protected:

float DensityToMass(float density) override;

public :

	explicit GrapheneFluid2D(SetUpParameters &input_parameters);
		~GrapheneFluid2D();

		/*!
		 * @brief Calculates @f$\Delta x@f$ and imposes Courant–Friedrichs–Lewy condition to @f$\Delta t@f$
		 *
		 * The method retrieves the spatial discretization @f$\Delta x = L / N_x @f$ and then sets the time step as @f$\Delta t = \Delta x/\lambda  @f$
		 * where @f$\lambda@f$ is given by
		 * @f{*eqnarray}
		   \lambda =1.2v_F \quad&\mathrm{if}\quad S<0.36v_F\\
		   \lambda =1.97S + 0.5 v_F \quad&\mathrm{otherwise} @f}
		 */
		void CflCondition() override;
		/*!
		 * @brief Sets the total simulation time in units of @f$v_0/L@f$
		 *
		 * The method predicts and set the maximum time for a typical Dyakonov-Shur instability simulation such that:
		  -# It is larger than the usual saturation time
		  -# After the saturation at least 10 oscillations occur

		  Our tests concluded heuristically that such criteria can be met setting @f$T_{max}=5+0.02S+20/S@f$
		 */
		void SetSimulationTime() override;

		/*!
		 * @brief Converts the mass density flux to velocity on the entire simulation grid.
		 *
		 * Since the mass of the fluid element is not a constant the in the graphene electronic fluid, one needs to perform the transformation
		   @f[ \vec{v} = \frac{\vec{p}}{n^{3/2}} @f]
		 * */
		void MassFluxToVelocity() override; // Converts the mass density flux back to velocity, in graphene  v = p n^{-3/2}
		/*Override fluxes and sources to specifics of graphene physics*/


//		float DensityFluxX(float n, float flx_x, float flx_y,float mass, float s) override;    ///< density equation (continuity equation) conserved flux X component
//		float DensityFluxY(float n, float flx_x, float flx_y,float mass, float s) override;    ///< density equation (continuity equation) conserved flux Y component
		float DensitySource(float n, float flx_x, float flx_y, float mass, float s)override;   ///< density equation (continuity equation) source term
//		float XMomentumFluxX(float n, float flx_x, float flx_y, float mass, float s) override;  ///< velocity X component equation (momentum equation) conserved flux X component
//		float XMomentumFluxY(float n, float flx_x, float flx_y, float mass, float s) override;  ///< velocity X component equation (momentum equation) conserved flux Y component
		float XMomentumSource(float n, float flx_x, float flx_y, float mass, float s)override; ///< velocity X component equation (momentum equation) source term
//		float YMomentumFluxX(float n, float flx_x, float flx_y, float mass, float s) override;  ///< velocity Y component equation (momentum equation) conserved flux X component
//		float YMomentumFluxY(float n, float flx_x, float flx_y, float mass, float s) override;  ///< velocity Y component equation (momentum equation) conserved flux Y component
		float YMomentumSource(float n, float flx_x, float flx_y, float mass, float s)override; ///< velocity y component equation (momentum equation) source term


	float DensityFluxX(GridPoint p, char side ) override; ///< density equation (continuity equation) conserved flux X component
	float DensityFluxY(GridPoint p, char side ) override; ///< density equation (continuity equation) conserved1 flux Y component

	float XMomentumFluxX(GridPoint p, char side ) override; ///< velocity X component equation (momentum equation) conserved flux X component
	float XMomentumFluxY(GridPoint p, char side ) override; ///< velocity X component equation (momentum equation) conserved flux Y component

	float YMomentumFluxX(GridPoint p, char side ) override; ///< velocity Y component equation (momentum equation) conserved flux X component
	float YMomentumFluxY(GridPoint p, char side ) override; ///< velocity Y component equation (momentum equation) conserved flux Y component


	void MagneticSourceSemiAnalytic(); // Semi analytic method for the magnetic interaction
		//void MagneticSourceFtcs();  // Forward Time Centered Space method for the magnetic interaction
};

#include "Fluid2DLib.h"

#endif //GRAPHENEFLUID2DLIB_H
