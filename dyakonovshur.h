
#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif


#ifndef M_EULER
#    define M_EULER 2.71828182845905
#endif

// Time prefactor in pico seconds
#ifndef PRE_TIM
#    define PRE_TIM 3.3333333333333
#endif

// Power prefactor in pico Watts 
#ifndef PRE_POW
#    define PRE_POW 9.53395
#endif

// Kinetic Energy prefactor in pico Joules
#ifndef PRE_KIN
#    define PRE_KIN 0.0168226
#endif

// Current prefactor in mili Amperes
#ifndef PRE_CUR
#    define PRE_CUR 0.480653
#endif

// VDS prefactor in Volts
#ifndef PRE_VDS
#    define PRE_VDS 18.0951
#endif


#ifndef Cspeed
#    define Cspeed 1000.0
#endif




float RealFreq(float s, float v, float L, int n);
float ImagFreq(float s, float v, float L);
void BoundaryCond(int type, int N, float * n, float * v);
void InitialCondSine(int N, float dx,  float * n, float * p);
void InitialCondRand(int N, float dx,  float * n, float * p);
float DtElectricDipole(int N,float dx, float * J);
float TotalCurrent(int N,float dx, float * n,float * v);
float TotalElectricDipole(int N,float dx, float * n);
float KineticEnergy(int N,float dx, float * n, float * v);
float RMS(int N, float dt, float * f);
float AVG(int N, float dt, float * f);
float GaussKernel(int n , float t);
float D_GaussKernel(int n , float t);
void ConvolveGauss(int type, float M, float t, float * in, float * out, int size);
void AverageFilter(float * vec_in, float * vec_out, int size , int width );
void ExtremaFinding(float * vec_in, int N, float S, float dt,float & sat, float  & tau, float & error, std::string extremafile);
void Autocorrelation(float * out_gamma ,float * in , int crop, int size);
void ShockFinding(float * in, int N, float t , float dx,  std::string shockfile);
void TimeDerivative(int size_rows,int size_cols, float dt,float ** f_in , float ** df_out );
void SpaceDerivative(int size_rows,int size_cols, float dt,float ** f_in , float ** df_out );
	 


float R(float x , float y , float z, float X , float Y , float Z );
float Retarded_time(float time, float x , float y , float z, float X , float Y , float Z );	 
void JefimenkoEMfield(int XDIM, int YDIM, float dx, float dy, float dt, float Xpos, float Ypos, float Zpos,  float ** rho, float ** rho_dot, float ** cur, float ** cur_dot, float Time , float  * E_out , float  * B_out, float  * S_out   );
