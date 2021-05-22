//
// Created by pcosme on 23/12/2020.
//

#include "DyakonovShurBoundaryLib.h"
#include "includes/BoundaryLib.h"
#include "RobinBoundaryLib.h"

void RobinBoundaryCondition::SlipLength(Fluid2D &fluid_class, float slip_length) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	float dy=fluid_class.GetDy();
	for (int i=0; i < nx; i++){
		float aux_top,aux_bot, dn_top,dn_bot,l_top,l_bot;
		int bottom=i; //i+0*nx
		int top= i + (ny - 1) * nx;
		dn_bot = (-3.0f*fluid_class.Den[i+0*nx]+4.0f*fluid_class.Den[i+1*nx]-1.0f*fluid_class.Den[i+2*nx])/(2.0f*dy);
		dn_top = (-3.0f*fluid_class.Den[i+(ny - 1)*nx]+4.0f*fluid_class.Den[i+(ny - 2)*nx]-1.0f*fluid_class.Den[i+(ny - 3)*nx])/(2.0f*dy);
		dn_bot = dn_bot/sqrt(fluid_class.Den[bottom]);
		dn_top = dn_top/sqrt(fluid_class.Den[top]);
		l_top = slip_length/(1.0f+slip_length*dn_top);
		l_bot = slip_length/(1.0f-slip_length*dn_bot);
		aux_top = l_top/(2.0f*dy+3.0f*l_top);
		aux_bot = l_bot/(2.0f*dy+3.0f*l_bot);
		fluid_class.Den[bottom] = fluid_class.Den[bottom + nx];
		fluid_class.Den[top] = fluid_class.Den[top - nx];
		fluid_class.FlxY[top] =  0.0f;
		fluid_class.FlxY[bottom] = 0.0f;
		fluid_class.FlxX[bottom] = aux_bot*(4.0f*fluid_class.FlxX[i+1*nx]-1.0f*fluid_class.FlxX[i+2*nx]);
		fluid_class.FlxX[top] = aux_top*(4.0f*fluid_class.FlxX[i+(ny-2)*nx]-1.0f*fluid_class.FlxX[i+(ny-3)*nx]);
	}
}

void RobinBoundaryCondition::SlipLengthTop(Fluid2D &fluid_class, float slip_length) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	float dy=fluid_class.GetDy();
	for (int i=0; i < nx; i++){
		float aux_top, dn_top,l_top;
		int top= i + (ny - 1) * nx;
		dn_top = (-3.0f*fluid_class.Den[i+(ny - 1)*nx]+4.0f*fluid_class.Den[i+(ny - 2)*nx]-1.0f*fluid_class.Den[i+(ny - 3)*nx])/(2.0f*dy);
		dn_top = dn_top/sqrt(fluid_class.Den[top]);
		l_top = slip_length/(1.0f+slip_length*dn_top);
		aux_top = l_top/(2.0f*dy+3.0f*l_top);
		fluid_class.Den[top] = fluid_class.Den[top - nx];
		fluid_class.FlxY[top] =  0.0f;
		fluid_class.FlxX[top] = aux_top*(4.0f*fluid_class.FlxX[i+(ny-2)*nx]-1.0f*fluid_class.FlxX[i+(ny-3)*nx]);
	}
}

void RobinBoundaryCondition::SlipLengthBottom(Fluid2D &fluid_class, float slip_length) {
	int nx=fluid_class.SizeX();
	float dy=fluid_class.GetDy();
	for (int i=0; i < nx; i++){
		float aux_bot, dn_bot, l_bot;
		//int bottom=i; //i+0*nx
		dn_bot = (-3.0f*fluid_class.Den[i+0*nx]+4.0f*fluid_class.Den[i+1*nx]-1.0f*fluid_class.Den[i+2*nx])/(2.0f*dy);
		dn_bot = dn_bot/sqrt(fluid_class.Den[i]);
		l_bot = slip_length/(1.0f-slip_length*dn_bot);
		aux_bot = l_bot/(2.0f*dy+3.0f*l_bot);
		fluid_class.Den[i] = fluid_class.Den[i + nx];
		fluid_class.FlxY[i] = 0.0f;
		fluid_class.FlxX[i] = aux_bot*(4.0f*fluid_class.FlxX[i+1*nx]-1.0f*fluid_class.FlxX[i+2*nx]);
	}
}