//
// Created by pcosme on 23/12/2020.
//

#ifndef DIRICHELETBOUNDARYLIB_H
#define DIRICHELETBOUNDARYLIB_H

#include <H5Cpp.h>
#include "TethysBaseLib.h"
#include "TethysMathLib.h"
#include "Fluid1DLib.h"
#include "Fluid2DLib.h"
#include "BoundaryLib.h"
#include "DiricheletBoundaryLib.h"



/*!
 * @brief Class for constant Dirichelet boundary conditions
 *
 * The class implements boundary conditions where the variables of the system are set to a constant value along the boundary. For the fluid this can translate in a variety of situations and physical scenarios.
 *
 * */
class  DirichletBoundaryCondition : public BoundaryCondition
{
	public:
	static void Density(Fluid1D& fluid_class, float left, float right);                                    ///< Fixed density at boundary n(x=0)=left and n(x=L)=right
	static void Density(Fluid2D& fluid_class, float left, float right, float top, float bottom);           ///< Fixed density at boundary n(x=0)=left, n(x=L)=right, n(y=0)=bottom, n(y=W)=top
	static void VelocityX(Fluid1D& fluid_class, float left, float right);                                  ///< Fixed Velocity at boundary V(x=0)=left and V(x=L)=right
	static void MassFluxX(Fluid2D& fluid_class, float left, float right, float top, float bottom);         ///< Fixed mass density flux x component at boundary px(x=0)=left, px(x=L)=right, px(y=0)=bottom, px(y=W)=top
	static void MassFluxY(Fluid2D& fluid_class, float left, float right, float top, float bottom);         ///< Fixed mass density flux y component at boundary py(x=0)=left, py(x=L)=right, py(y=0)=bottom, py(y=W)=top
	static void Jet(Fluid2D& fluid_class, float left, float left_width, float right, float right_width);   ///< Jet configuration i.e. fixed flux x component at a portion of given with around the center of the edges x=0 and x=L. Useful to study turbulence onset
	static void DensityRight(Fluid2D& fluid_class, float right);       ///< Fixed density at boundary n(x=L)=right
	static void MassFluxXRight(Fluid2D& fluid_class, float right);     ///< Fixed mass density flux X component at boundary px(x=L)=right
	static void MassFluxYRight(Fluid2D& fluid_class, float right);     ///< Fixed mass density flux Y component at boundary py(x=L)=right
	static void DensityLeft(Fluid2D& fluid_class, float left);         ///< Fixed density at boundary n(x=0)=left
	static void MassFluxXLeft(Fluid2D& fluid_class, float left);       ///< Fixed mass density flux X component at boundary px(x=0)=left
	static void MassFluxYLeft(Fluid2D& fluid_class, float left);       ///< Fixed mass density flux Y component at boundary py(x=0)=left
	static void DensityTop(Fluid2D& fluid_class, float top);           ///< Fixed density at boundary n(y=W)=top
	static void MassFluxXTop(Fluid2D& fluid_class, float top);         ///< Fixed mass density flux X component at boundary px(y=W)=top
	static void MassFluxYTop(Fluid2D& fluid_class, float top);         ///< Fixed mass density flux Y component at boundary py(y=W)=top
	static void DensityBottom(Fluid2D& fluid_class, float bottom);     ///< Fixed density at boundary n(y=0)=bottom
	static void MassFluxXBottom(Fluid2D& fluid_class, float bottom);   ///< Fixed mass density flux x component at boundary px(y=0)=bottom
	static void MassFluxYBottom(Fluid2D& fluid_class, float bottom);   ///< Fixed mass density flux Y component at boundary py(y=0)=bottom
    static void Temperature(Fluid2D& fluid_class, float left, float right, float top, float bottom);           ///< Fixed density at boundary n(x=0)=left, n(x=L)=right, n(y=0)=bottom, n(y=W)=top
    static void TemperatureRight(Fluid2D& fluid_class, float right);       ///< Fixed density at boundary n(x=L)=right
    static void TemperatureLeft(Fluid2D& fluid_class, float left);         ///< Fixed density at boundary n(x=0)=left
    static void TemperatureTop(Fluid2D& fluid_class, float top);           ///< Fixed density at boundary n(y=W)=top
    static void TemperatureBottom(Fluid2D& fluid_class, float bottom);     ///< Fixed density at boundary n(y=0)=botto
};


#endif //DIRICHELETBOUNDARYLIB_H
