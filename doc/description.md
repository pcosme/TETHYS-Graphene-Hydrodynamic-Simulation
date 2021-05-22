@page  description Software Description
@tableofcontents


@section purpose Purpose and scope

The main objective of this code is to simulate and analyse a hydrodynamic description of electronic conduction on (mono-layer) graphene. At the present time limiting the simulations (we refer to one spatial dimension simulation when the geometry and symmetry of the problem let one foreseen that the dynamics of the system will develop solely along the direction of the main electronic current i.e. along the GFET channel and that the quantities and derivatives along the transverse direction of the plane can be discarded accordingly. If, on the contrary, the fluid is expected to display a non parallel flow or transverse variations a complete 2D simulation is required.), both 1D+1 and 2D+1, to the case of gated graphene system. The physical description of such problem can be found at \cite Cosme2019 . In the bidimensional simulation is possible to add the effects of shear viscosity and/or uniform and static magnetic field along the transverse direction, while the implementation of odd viscosity \cite Avron1998OddViscosity , also commonly dubbed as *Hall viscosity*, is ongoing.

As a typical hydrodynamic simulation it solves a conservation law equation to obtain the density, @f$n(x,y,t)@f$, and velocity, @f$\vec{v}(x,y,t)@f$, profiles. The equations governing this problem can be adimensionalised redefining the variables in terms of the channel length @f$L@f$ and characteristic velocity and density @f$v_0@f$ and @f$n_0@f$, that can be taken as the equilibrium or steady state values, writing
@f[
x^\ast\equiv x/L\quad t^\ast\equiv tv_0/L\quad \omega^\ast\equiv \omega L/v_0\quad v^\ast\equiv v/v_0\quad n^\ast\equiv n/n_0 
@f]
thus, from this point forward all the quantities are implied to be dimensionless even though the asterik sign  will be dropped for simplicity.

The complete model that the code aims to simulated is given by the equations


@f[
\frac{\partial n}{\partial t} +\bm{\nabla}\!\!\cdot\! n\mathbf{v} = 0   
@f]
and
@f[
\frac{\partial \vec{v}}{\partial t} + \frac{1}{2}(\vec{v}\cdot \bm{\nabla} )\vec{v}  +\frac{v_F^2}{2n}\bm{\nabla}n+\frac{S^2}{\sqrt{n}}\bm{\nabla}n-\\-\nu_s\nabla^2\vec{v}-\nu_o\nabla^2\vec{v}^\dagger+\omega_c\vec{v}\times\vec{\hat{z}}=\frac{1-\vec{v}}{\tau},\label{eq:eulerequations2}    
@f]

where @f$\vec{v}^\dagger=[-v_y,v_x]^\mathsf{T}@f$ refers to the rotated velocity vector and we consider a magnetic field @f$\vec{B}=B_0\vec{\hat{z}}@f$.  All the parameters are detailed at Table ?. Evidently, some aspects of the system, as the odd viscosity and the magnetic terms, can only have meaning in a 2D simulation.
The parameter @f$S@f$ has units of a velocity and can be interpreted as a sound speed of the electron fluid (not to be confused with the intrinsic sound speed of phonons in graphene). Moreover, the ratio @f$S/v_0@f$ will play a crucial role determining the properties of the hydrodynamic system, similarly to the Froude number in fluid dynamics. For typical parameters of a graphene FET, @f$S/v_0@f$ scales up to a few tens.
While the system can easily be casted in a conservation form form the one-dimensional case, the two-dimensional scenario poses some technical difficulties in writing as a flux law. For that reason, we resorted to the mass flux density @f$\vec{p}=m^\star n\vec{v}@f$ which allows us to rewrite the model in the equivalent form   
@f[
\frac{\partial n}{\partial t}+\bm{\nabla}\!\cdot\!\frac{\vec{p}}{\sqrt{n}}=0 \text{  and}    
\label{eq:moment_1}
@f]
@f[
\frac{\partial \vec{p}}{\partial t}+\bm{\nabla}\!\cdot\!\left( \frac{ \vec{p}\otimes \vec{p}}{{n}^{3/2}} + v_F^2\frac{n^{3/2}}{3}\mathds{1}+S^2\frac{{n}^2}{2}\mathds{1}\right)+\\
-\nu_s\nabla^2\frac{\vec{p}}{n^{3/2}} -\nu_o\nabla^2\frac{\vec{p}^\dagger}{n^{3/2}}
+\omega_c\frac{\vec{p}\times\vec{\hat{z}}}{\sqrt{n}}=0. \label{eq:moment_2}
@f]

that was implemented in the 2D algorithms.
 

| Parameter   |      Description    |   Relations  | Requirements | Typical values | 
| :---:   |      :---      |  --- | :---:| ---  |
| @f$ v_F @f$   |  Fermi velocity @f${}^{(a)}@f$   |  |   |  @f$1\times 10^6 ms^{-1}@f$  |
| @f$E_F@f$   |  Fermi level    | @f$= \hbar v_F\sqrt{\pi n}@f$  |  @f$ k_BT\ll E_F \ll t \simeq 3 eV@f$ |  @f$0.0658 \sqrt{n} eV \dagger@f$ |
| @f$\mu@f$  |   Carrier mobility  |   |   | @f$2000~5000 cm^2/Vs@f$   |
| @f$\tau@f$  |   Relaxation time  | @f$= \lambda_i/v_F@f$  |   | @f$\sim 5ps   (0.1 \sim 10ps)@f$ |
| @f$1/\tau@f$  |   Collision frequency@f${}^{(a)}@f$  |   |   | @f$\sim0.2THz @f$ |
| @f$\lambda_i@f$  |  Mean free path (impurities)     | |  |  @f$\sim5\mu m@f$ |
| @f$\lambda_{ee}@f$  | Mean free path (electrons)       ||   | @f$\sim10nm@f$  |
| @f$ v_{sat} @f$  |   Saturation velocity  |     @f$\sim0.5v_F@f$| | @f$0.5\times10^6 ms^{-1}@f$   |
| @f$\nu_s@f$  |   Shear kinematic viscosity@f${}^{(a)}@f$     ||   |  @f$\sim 0.2  m^2s^{-1}@f$  |
| @f$\nu_o@f$  |   Odd kinematic viscosity     || broken time symmetry (@f$B\neq0@f$)  | @f$ \lesssim 0.1  m^2s^{-1}@f$  |
| @f$C_g@f$ | Capacitance per area   |  @f$= \varepsilon/d_0@f$|  | @f$0.1~0.6 \muFcm^{-2}@f$  |
| @f$L@f$ |   Channel length@f${}^{(b)}@f$  |  | @f$\lambda_{ee} < L <  \lambda_i @f$   | |
| @f$W@f$ |  Channel width   |  | @f$\lambda_{ee} < W < \lambda_i  @f$  | |
| @f$x:y@f$ |  Aspect ratio@f${}^{(a)}@f$   | = @f$L/W@f$ |    | |
| @f$v_0@f$ |   Mean drift velocity@f${}^{(b)}@f$  |  |  @f$\ll v_F and < v_sat@f$| @f$~1x10^5 ms^{-1}@f$|
| @f$n_0@f$ |   Mean carrier density@f${}^{(b)}@f$  |  | | @f$~ 1\times10^{12} cm^{-2}@f$ |
| @f$m^\star@f$ | Effective mass  | @f$= \hbar\sqrt{\pi n}/v_F@f$ | | @f$\sim0.02m_e@f$ |
| @f$S@f$ |  Plasmons velocity@f${}^{(a)}@f$   | @f$= ev_F n_0^{1/4}/\hbarC_g^{1/2}\pi^{1/2}@f$    |  @f$> v_F@f$ |@f$1.17\times10^6 n_0^{1/4}C_g^{-1/2}  ms^{-1}\dagger @f$|
| @f$Re@f$  |    Reynolds number    | @f$= v_0L/\nu_s@f$   |   | |
| @f$t_0@f$ |  Time scale@f${}^{(b)}@f$  | = @f$L/v_0@f$   |   | |
| @f$\omega_0@f$ |  Frequency scale@f${}^{(b)}@f$ |@f$ = v_0/L@f$   |   |  |
| @f$\omega_c@f$ |  Cyclotron frequency@f${}^{(a)}@f$  |@f$ =eB/m^\star@f$   |  @f$B\ll1 T @f$| @f$\ll 9 THz@f$ |
| @f$\omega_{DS}@f$ |   Dyakonov-Shur frequency  |  @f$\sim \piS/2L@f$  |   |  |
| @f$\gamma_{DS}@f$ |  Dyakonov-Shur growth rate   | @f$\sim0.75 v_0/L@f$ | > @f$1/\tau@f$  |  |
| @f$U_0@f$ |   Gate mean voltage  | @f$= en_0/C_g@f$  | | @f$160 n_0/C_g mV\dagger@f$ |
| @f$I_0@f$ |   Typical channel current  | @f$= en_0v_0W@f$ | @f$<1.3 mA/cm for L=25 \mu m@f$  | @f$160 n_0v_0W \muA \ddag@f$ |

<p>@f${}^{(a)}@f$ User defined parameter</p>
<p>@f${}^{(b)}@f$ Normalisation parameter</p>
<p>@f$\dagger@f$ If the density  @f$n_0@f$ is given in @f$10^12 cm^{-2}@f$  and the capacitance @f$C_g@f$ in @f$\mu Fcm^{-2}@f$ </p>
<p>@f$\ddag@f$ If the density  @f$n_0 is@f$ given in @f$10^12 cm^{-2}@f$, the drift velocity @f$v_0@f$ in @f$10^5 ms^{-1}@f$ and the channel width W in @f$\mu m@f$ </p>


@subpage numericalmethods
