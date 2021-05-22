//
// Created by pcosme on 23/12/2020.
//


#ifndef ROBINBOUNDARYLIB_H
#define ROBINBOUNDARYLIB_H
#include <H5Cpp.h>
#include "TethysBaseLib.h"
#include "TethysMathLib.h"
#include "Fluid1DLib.h"
#include "Fluid2DLib.h"
#include "DiricheletBoundaryLib.h"

/*!
 * @brief Class for Robin type boundary conditions.
 *
 * This type of boundary conditions are characterized by a linear combination of the value and the normal derivative of the quantities at the boundary, i.e.
 * @f[ a U +  b \frac{\partial U}{\partial \xi}= c\quad {\rm at}\quad \partial \Omega\quad {\rm with}\quad a,b,c\in\mathbb{R}  @f]
 *
 * In the context of the fluid simulation it is useful for setting the _slip length_ conditions at the interfaces, where the fluid velocity is set as
 *
 * @f[ v_x = \pm\ell \frac{\partial v_x}{\partial y} \quad v_y=0@f]
 * */
class  RobinBoundaryCondition : public DirichletBoundaryCondition
{
public:
	static void SlipLength(Fluid2D& fluid_class,float slip_length);        ///< slip length boundary condition Vx = +-l dVx/dy and Vy=0 at y=0 and y=W
	static void SlipLengthTop(Fluid2D& fluid_class,float slip_length);     ///< slip length boundary condition Vx = -l dVx/dy and Vy=0 at y=W
	static void SlipLengthBottom(Fluid2D& fluid_class,float slip_length);  ///< slip length boundary condition Vx = l dVx/dy and Vy=0 at y=0
};

#endif //ROBINBOUNDARYLIB_H
