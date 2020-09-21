#ifndef BOUNDARYLIB_H
#define BOUNDARYLIB_H

#include <H5Cpp.h>
#include "TethysLib.h"
#include "Tethys1DLib.h"
#include "Tethys2DLib.h"

class BoundaryCondition {
	public :
	static void XFree(Fluid1D& fluid_class);
	void XFree(Fluid2D& fluid_class);
	void XFreeLeft(Fluid2D& fluid_class);
	void XFreeRight(Fluid2D& fluid_class);
	void XPeriodic(Fluid1D& fluid_class);
	void XPeriodic(Fluid2D& fluid_class);
	void YFree(Fluid2D& fluid_class);
	void YFreeTop(Fluid2D& fluid_class);
	void YFreeBottom(Fluid2D& fluid_class);
	void YPeriodic(Fluid2D& fluid_class);
	void YClosedFreeSlip(Fluid2D& fluid_class);
	void YClosedNoSlip(Fluid2D& fluid_class);

	//class DyakonovShur;
	//class Dirichlet;
};	



class  DirichletBoundaryCondition : public BoundaryCondition
{
	public: 	
	void Density(Fluid1D& fluid_class, float left, float right);
	void Density(Fluid2D& fluid_class, float left, float right, float top, float bottom);
	void VelocityX(Fluid1D& fluid_class, float left, float right);
	void MassFluxX(Fluid2D& fluid_class, float left, float right, float top, float bottom);
	void MassFluxY(Fluid2D& fluid_class, float left, float right, float top, float bottom);
	void Jet(Fluid2D& fluid_class, float left, float left_width, float right, float right_width);


	void DensityRigth(Fluid2D& fluid_class, float right);
	void MassFluxXRight(Fluid2D& fluid_class, float right);
	void MassFluxYRight(Fluid2D& fluid_class, float right);

	void DensityLeft(Fluid2D& fluid_class, float left);
	void MassFluxXLeft(Fluid2D& fluid_class, float left);
	void MassFluxYLeft(Fluid2D& fluid_class, float left);

	void DensityTop(Fluid2D& fluid_class, float top);
	void MassFluxXTop(Fluid2D& fluid_class, float top);
	void MassFluxYTop(Fluid2D& fluid_class, float top);

	void DensityBottom(Fluid2D& fluid_class, float bottom);
	void MassFluxXBottom(Fluid2D& fluid_class, float bottom);
	void MassFluxYBottom(Fluid2D& fluid_class, float bottom);
};

class  DyakonovShurBoundaryCondition : public DirichletBoundaryCondition
{
public:
	void DyakonovShurBC(GrapheneFluid1D& fluid_class);
	void DyakonovShurBC(GrapheneFluid2D& fluid_class);
};


#endif

