/************************************************************************************************\
* 2022 Pedro Cosme , João Santos , Ivan Figueiredo and Diogo Simões                              *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/

#include "includes/BoundaryLib.h"
#include "includes/FeedbackBoundaryLib.h"

FeedbackBoundaryCondition::FeedbackBoundaryCondition(float time_delay, float delta){
    count=0;
    Nsteps=(int)(time_delay/delta);
    Curr = new float[Nsteps];
    Dens = new float[Nsteps];
    for(int i = 0; i < Nsteps; ++i){
        Curr[i]=0;
        Dens[i]=0;
    }
}

FeedbackBoundaryCondition::~FeedbackBoundaryCondition()=default;

void FeedbackBoundaryCondition::VoltageFeedbackBc(GrapheneFluid1D& fluid_class, const float* Trans, float intens, float omega, float t) {
	int nx=fluid_class.SizeX();
    float Vi = fluid_class.Umain[nx-1].n();
    float Ii = fluid_class.Umain[nx-1].v()*fluid_class.Umain[nx-1].n();
    
    float Vf = Trans[0]*(Vi-1)+Trans[1]*(Ii-1);
    float If = Trans[2]*(Vi-1)+Trans[3]*(Ii-1);
    
    fluid_class.Umain[0].n()+=Vf+intens*cos(omega*t);
    fluid_class.Umain[0].v()+=If/fluid_class.Umain[0].n();
}

void FeedbackBoundaryCondition::CurrentFeedbackBc(GrapheneFluid1D& fluid_class, const float* Trans, float intens, float omega, float t) {
	int nx=fluid_class.SizeX();
    float Vi = fluid_class.Umain[nx-1].n();
    float Ii = fluid_class.Umain[nx-1].v()*fluid_class.Umain[nx-1].n();
    
    float Vf = Trans[0]*(Vi-1)+Trans[1]*(Ii-1);
    float If = Trans[2]*(Vi-1)+Trans[3]*(Ii-1);
    
    fluid_class.Umain[0].n()+=Vf;
    fluid_class.Umain[0].v()+=(If+intens*cos(omega*t))/fluid_class.Umain[0].n();
}

void FeedbackBoundaryCondition::VoltageDelayFeedbackBc(GrapheneFluid1D& fluid_class, const float* Trans, float intens, float omega, float t) {
	int nx=fluid_class.SizeX();
    Dens[count%Nsteps] = fluid_class.Umain[nx-1].n();
    Curr[count%Nsteps] = fluid_class.Umain[nx-1].v()*fluid_class.Umain[nx-1].n();
    count=(count+1)%Nsteps;
    
    float Vf = Trans[0]*(Dens[count]-1)+Trans[1]*(Curr[count]-1);
    float If = Trans[2]*(Dens[count]-1)+Trans[3]*(Curr[count]-1);

	fluid_class.Umain[0].n()+=Vf+intens*cos(omega*t);
	fluid_class.Umain[0].v()+=If/fluid_class.Umain[0].n();
}

void FeedbackBoundaryCondition::CurrentDelayFeedbackBc(GrapheneFluid1D& fluid_class, const float* Trans, float intens, float omega, float t) {
	int nx=fluid_class.SizeX();
    Dens[count%Nsteps] = fluid_class.Umain[nx-1].n();
    Curr[count%Nsteps] = fluid_class.Umain[nx-1].v()*fluid_class.Umain[nx-1].n();
    count=(count+1)%Nsteps;

    float Vf = Trans[0]*(Dens[count]-1)+Trans[1]*(Curr[count]-1);
    float If = Trans[2]*(Dens[count]-1)+Trans[3]*(Curr[count]-1);

	fluid_class.Umain[0].n()+=Vf;
	fluid_class.Umain[0].v()+=(If+intens*cos(omega*t))/fluid_class.Umain[0].n();
}