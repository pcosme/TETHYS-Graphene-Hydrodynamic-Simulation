#ifndef BOUNDARYLIB_H
#define BOUNDARYLIB_H

#include <H5Cpp.h>
#include "TethysLib.h"
#include "TethysMathLib.h"
#include "Tethys1DLib.h"
#include "Tethys2DLib.h"

/*Base class for general boundary conditions*/
class BoundaryCondition {
	public :
		void XFree(Fluid1D& fluid_class);           // open boundaries at x=0 and x=L for all variables and zero tangent velocity
		void XFree(Fluid2D& fluid_class);           // open boundaries at x=0 and x=L for all variables and zero tangent velocity
		void XFreeLeft(Fluid2D& fluid_class);       // open boundaries at x=0 for all variables and zero tangent velocity Vy=0
		void XFreeRight(Fluid2D& fluid_class);      // open boundaries at x=L for all variables and zero tangent velocity Vy=0
		void XPeriodic(Fluid1D& fluid_class);       // periodic boundaries u(x=0)=u(x=L) for all variables  and zero tangent velocity
		void XPeriodic(Fluid2D& fluid_class);       // periodic boundaries u(x=0)=u(x=L) for all variables and zero tangent velocity
		void YFree(Fluid2D& fluid_class);           // open boundaries at y=0 and y=W for all variables
		void YFreeTop(Fluid2D& fluid_class);        // open boundaries at y=W for all variables
		void YFreeBottom(Fluid2D& fluid_class);     // open boundaries at y=0 for all variables
		void YPeriodic(Fluid2D& fluid_class);       // periodic boundaries u(y=0)=u(y=W) for all variables
		void YClosedFreeSlip(Fluid2D& fluid_class); // zero flux across y=0 and y=W and free tangent velocity Vx
		void YClosedNoSlip(Fluid2D& fluid_class);   // zero flux across y=0 and y=W and zero tangent velocity Vx=0
};


/*Class for Dirichlet boundary conditions in the form U=const at the boundary*/
class  DirichletBoundaryCondition : public BoundaryCondition
{
	public: 	
	void Density(Fluid1D& fluid_class, float left, float right);                                    // Fixed density at boundary n(x=0)=left and n(x=L)=right
	void Density(Fluid2D& fluid_class, float left, float right, float top, float bottom);           // Fixed density at boundary n(x=0)=left, n(x=L)=right, n(y=0)=bottom, n(y=W)=top
	void VelocityX(Fluid1D& fluid_class, float left, float right);                                  // Fixed Velocity at boundary V(x=0)=left and V(x=L)=right
	void MassFluxX(Fluid2D& fluid_class, float left, float right, float top, float bottom);         // Fixed mass density flux x component at boundary px(x=0)=left, px(x=L)=right, px(y=0)=bottom, px(y=W)=top
	void MassFluxY(Fluid2D& fluid_class, float left, float right, float top, float bottom);         // Fixed mass density flux y component at boundary py(x=0)=left, py(x=L)=right, py(y=0)=bottom, py(y=W)=top
	void Jet(Fluid2D& fluid_class, float left, float left_width, float right, float right_width);   // Jet configuration i.e. fixed flux x component at a portion of given with around the center of the edges x=0 and x=L. Useful to study turbulence onset
	void DensityRight(Fluid2D& fluid_class, float right);       // Fixed density at boundary n(x=L)=right
	void MassFluxXRight(Fluid2D& fluid_class, float right);     // Fixed mass density flux X component at boundary px(x=L)=right
	void MassFluxYRight(Fluid2D& fluid_class, float right);     // Fixed mass density flux Y component at boundary py(x=L)=right
	void DensityLeft(Fluid2D& fluid_class, float left);         // Fixed density at boundary n(x=0)=left
	void MassFluxXLeft(Fluid2D& fluid_class, float left);       // Fixed mass density flux X component at boundary px(x=0)=left
	void MassFluxYLeft(Fluid2D& fluid_class, float left);       // Fixed mass density flux Y component at boundary py(x=0)=left
	void DensityTop(Fluid2D& fluid_class, float top);           // Fixed density at boundary n(y=W)=top
	void MassFluxXTop(Fluid2D& fluid_class, float top);         // Fixed mass density flux X component at boundary px(y=W)=top
	void MassFluxYTop(Fluid2D& fluid_class, float top);         // Fixed mass density flux Y component at boundary py(y=W)=top
	void DensityBottom(Fluid2D& fluid_class, float bottom);     // Fixed density at boundary n(y=0)=bottom
	void MassFluxXBottom(Fluid2D& fluid_class, float bottom);   // Fixed mass density flux x component at boundary px(y=0)=bottom
	void MassFluxYBottom(Fluid2D& fluid_class, float bottom);   // Fixed mass density flux Y component at boundary py(y=0)=bottom
};

/*Class for the asymmetric Dyakonov-Shur boundary conditions*/
class  DyakonovShurBoundaryCondition : public DirichletBoundaryCondition
{
public:
	void DyakonovShurBc(GrapheneFluid1D& fluid_class);  // Dyakonov-Shur boundary conditions 1D n(0)=1 n(L)V(L)=1
	void DyakonovShurBc(GrapheneFluid2D& fluid_class);  // Dyakonov-Shur boundary conditions 2D n(x=0)=1 n(x=L)Vx(x=L)=1 Vy(x=0)=0 Vy(L=0)=0
};

/*Class for Robin boundary conditions in the form aU+bU'=const at the boundary*/
class  RobinBoundaryCondition : public DirichletBoundaryCondition
{
public:
	void SlipLength(Fluid2D& fluid_class,float slip_length);        // slip length boundary condition Vx = l dVx/dy and Vy=0 at y=0 and y=W
	void SlipLengthTop(Fluid2D& fluid_class,float slip_length);     // slip length boundary condition Vx = l dVx/dy and Vy=0 at y=W
	void SlipLengthBottom(Fluid2D& fluid_class,float slip_length);  // slip length boundary condition Vx = l dVx/dy and Vy=0 at y=0
};

#endif

