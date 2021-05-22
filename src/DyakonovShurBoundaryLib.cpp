//
// Created by pcosme on 23/12/2020.
//

#include "includes/BoundaryLib.h"
#include "DyakonovShurBoundaryLib.h"

void DyakonovShurBoundaryCondition::DyakonovShurBc(GrapheneFluid1D& fluid_class) {
	int nx=fluid_class.SizeX();
	fluid_class.Den[0] = 1.0f;
	fluid_class.Vel[0] = fluid_class.Vel[1];
	fluid_class.Den[nx - 1] = fluid_class.Den[nx - 2];
	fluid_class.Vel[nx - 1] = 1.0f / fluid_class.Den[nx - 1];
}

void DyakonovShurBoundaryCondition::DyakonovShurBc(GrapheneFluid2D& fluid_class) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();

	for(int j=0; j < ny; j++) {
		fluid_class.Den[0 + j * nx] = 1.0f;            //constant density at x=0
		fluid_class.FlxX[0 + j * nx] =
				fluid_class.FlxX[1 + j * nx] * pow(fluid_class.Den[1 + j * nx], -1.5f);            //free flux at x=0
		fluid_class.FlxY[0 + j * nx] = 0.0f;                    //flux only on x at x=0

		fluid_class.Den[nx - 1 + j * nx] = fluid_class.Den[nx - 2 + j * nx];            //free density at x=L
		fluid_class.FlxX[nx - 1 + j * nx] = sqrt(fluid_class.Den[nx - 1 + j * nx]);    //constant current at x=L (flux equals mass)
		fluid_class.FlxY[nx - 1 + j * nx] = 0.0f;                    //idem at x=L

	}
}