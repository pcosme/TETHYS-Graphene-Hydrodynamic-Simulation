/************************************************************************************************\
* Copyright (c) 2020 Pedro Cosme and João Santos                                                 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/

#ifndef BOUNDARYLIB_H
#define BOUNDARYLIB_H

#include <H5Cpp.h>
#include "TethysBaseLib.h"
#include "TethysMathLib.h"
#include "Fluid1DLib.h"
#include "Fluid2DLib.h"


//TODO review the entire boundary class

/*Base class for general boundary conditions*/

/*!
 * @brief Base class for boundary conditions
 *
 * The class implements the more general and simple boundary conditions such as free flow or periodic ones.
 *
 * */
class BoundaryCondition {
	protected:
		static int * BottomEdge;
		static int * TopEdge;
		static float  Slope;
	public :
		static void SetTopEdge(Fluid2D &fluid_class);
		static void SetBottomEdge(Fluid2D &fluid_class);
		static void SetSlope(float boundary_slope);
		static float GetSlope() ;
		static void XFree(Fluid1D& fluid_class);           ///< open boundaries at x=0 and x=L for all variables and zero tangent velocity
		static void XFree(Fluid2D& fluid_class);           ///< open boundaries at x=0 and x=L for all variables and zero tangent velocity
		static void XFree(Fluid2D &fluid_class, int x_limit);       ///< open boundaries at x=0 for all variables and zero tangent velocity Vy=0
	//	static void XFreeRight(Fluid2D& fluid_class);      ///< open boundaries at x=L for all variables and zero tangent velocity Vy=0
		static void XPeriodic(Fluid1D& fluid_class);       ///< periodic boundaries u(x=0)=u(x=L) for all variables  and zero tangent velocity
		static void XPeriodic(Fluid2D& fluid_class);       ///< periodic boundaries u(x=0)=u(x=L) for all variables and zero tangent velocity
		static void YFree(Fluid2D& fluid_class);           ///< open boundaries at y=0 and y=W for all variables
		static void YFree(Fluid2D &fluid_class, int y_limit);        ///< open boundaries at y=W for all variables
	//	static void YFreeBottom(Fluid2D& fluid_class);     ///< open boundaries at y=0 for all variables
		static void YPeriodic(Fluid2D& fluid_class);       ///< periodic boundaries u(y=0)=u(y=W) for all variables
		static void YClosedFreeSlip(Fluid2D& fluid_class); ///< zero flux across y=0 and y=W and free tangent velocity Vx
		static void YClosedNoSlip(Fluid2D& fluid_class);   ///< zero flux across y=0 and y=W and zero tangent velocity Vx=0
	static void XFreeLeft(Fluid2D &fluid_class);

	static void XFreeRight(Fluid2D &fluid_class);

	static void YFreeTop(Fluid2D &fluid_class);

	static void YFreeBottom(Fluid2D &fluid_class);
};

#endif

