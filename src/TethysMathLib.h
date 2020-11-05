#ifndef TETHYSMATHLIB_H
#define TETHYSMATHLIB_H

#include <H5Cpp.h>
#include "TethysLib.h"

using namespace std;

/* Functions to implement the spatial variation of the sound velocity S(x) in 1D or S(x,y) in 2D
 * corresponding to a variation of substrat permitivitty or even the description of a multi gated system.
 * */
float Sound_Velocity_Anisotropy(float x, float s);
float Sound_Velocity_Anisotropy(float x, float y, float s);

void Convolve_Gauss(unsigned int type, unsigned int m, float t, const float * in, float * out, unsigned long size);
//void Convolve_Gauss(int type, float m, float t, float * in, float * out, int size);
constexpr float Gauss_Kernel(int pos , float t); //
constexpr float Gauss_Kernel_Derivative(int pos , float t); //
float Signal_Average(int n, float dt, const float * f);
float Integral_1_D(int n, float ds, const float * f);
float Integral_2_D(int n, int m, float dx, float dy, const float * f);
/* Average moving filter for the smoothing of 1D simulation, suppressing the spurious oscillations inherent to the 2nd order solver*/
void Average_Filter(const float * vec_in, float * vec_out, int size , int width );


#endif //TETHYSMATHLIB_H
