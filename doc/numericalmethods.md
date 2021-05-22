@page numericalmethods Numerical methods
@tableofcontents

The simulations of the model described before rely on the division of the problem into two separate routines -- one for the hyperbolic part (convective and conserved terms) and a second one for the parabolic (viscous and source terms). Viewing the model as having the form @f$\partial_t\vec{u}=\mathscr{L}_\text{conservation}(\vec{u})+\mathscr{L}_\text{dissipation}(\vec{u})@f$, so that an operator splitting method could be used, uncoupling the two terms, for instance, for the mass flux equation one has:
@f[
\frac{\partial \vec{p}}{\partial t}+\overbrace{\bm{\nabla}\!\cdot\!\left( \frac{ \vec{p}\otimes \vec{p}}{{n}^{3/2}} + v_F^2\frac{n^{3/2}}{3}\mathds{1}+S^2\frac{{n}^2}{2}\mathds{1}\right)}^{\mathscr{L}_\text{conservation}(\vec{p})}+\\
-\underbrace{\nu_s\nabla^2\frac{\vec{p}}{n^{3/2}} -\nu_o\nabla^2\frac{\vec{p}^\dagger}{n^{3/2}}
+\omega_c\frac{\vec{p}\times\vec{\hat{z}}}{\sqrt{n}}}_{\mathscr{L}_\text{dissipation}(\vec{p})}=0.
@f]

Then, if the solution, after a time @f$t@f$, for the conservation and dissipation operators is well approximated by the numerical evolution operators @f$\Lambda^t_\text{con}@f$ and @f$\Lambda^t_\text{dis}@f$, respectively, the solution of the complete system is given by the Godunov splitting. This method is generally first order in time. One can, of course, invert the order of the operators and in fact the splitting error depends on the commutator @f$[\mathscr{L}_\text{dis},\mathscr{L}_\text{con}]@f$. For second order accuracy one can use Strang splitting where @f$\vec{u}^{t_0+t}\approx\Lambda^{t/2}_\text{dis}\circ\Lambda^t_\text{con}\circ\Lambda^{t/2}_\text{dis}(\vec{u}^{t_0})@f$. However, the performed test so far have not find significant accuracy or stability improvement when using the more computation time costly Strang splitting so we recommend to apply the simple Godunov algorithm:
@f[
\vec{u}^{t_0+t}\approx\Lambda^t_\text{dis}\circ\Lambda^t_\text{con}(\vec{u}^{t_0}).
@f]

The hyperbolic operator was solved with a finite volume Lax-Wendroff type method \cite Hirsch2007 \cite LeVeque1992  --  the two step Richtmyer scheme for nonlinear systems \cite LeVeque1992 .  Such scheme was then elected to implement, given it is computationally light, as it is not required to explicitly evaluate the system Jacobian at each step, yet capable of returning accurate simulations. It consists in two distinct steps, a predictor where the mid nodes are calculated and then a corrector step that updates the simulated quantities at the central nodes.

Whereas the dissipation and source terms were then computed with a forward time centred space (FTCS) scheme of finite differences. This choice of a lightweight parabolic solver has impact on the viable set of parameters for simulation, nonetheless, it allows to keep the maximum computing time for a complete 2D simulation bellow 150 min.

@section onedrich One-dimensional Richtmyer method

In the onedimensional case, the hyperbolic system of equations for density and velocity fields is written in a conservation form:
@f[
\frac{\partial}{\partial t}\begin{bmatrix}n\\v\end{bmatrix}+\frac{\partial}{\partial x}
\begin{bmatrix}      nv\\ \frac{v^2}{4}+\frac{v_F^2}{2}\log n + 2S^2\sqrt{n}\\\end{bmatrix} = 0 \iff \frac{\partial}{\partial t}\vec{u}+\frac{\partial}{\partial x}\mathbf{F(u)}=0, \label{eq:conserv}
@f]

and was implemented for @f$n@f$ an @f$v@f$ according to

@f{align}
\vec{u}_{i+\sfrac{1}{2}}^{k+\sfrac{1}{2}} &= \frac{1}{2}\left(\vec{u}_{i+1}^k + \vec{u}_{i}^k\right) - \frac{\Delta t}{2\Delta x}\left( \mathbf{F}_{i+1}^k - \mathbf{F}_{i}^k\right)\label{eq:richt_1}\\
\vec{u}_i^{k+1} &= \vec{u}_i^k - \frac{\Delta t}{\Delta x} \left( \mathbf{F}_{i+\sfrac{1}{2}}^{k+\sfrac{1}{2}} - \mathbf{F}_{i-\sfrac{1}{2}}^{k+\sfrac{1}{2}}\right)\label{eq:richt_2}
@f}
using @f$\vec{u}@f$ and @f$\mathbf{F}(\vec{u})@f$ as defined by equation and space (indicated at subscript indices) and time (indicated at superscript indices) discretisation where @f$t=k\Delta t@f$ and @f$t=i\Delta x@f$. Since such scheme is second order both in time and space it will not introduce spurious diffusion in the solution as a first order scheme would. It will, however, introduce artificial oscillations, contiguous to discontinuities, which can be corrected by the application of a moving average smoothing filter.

To comply with Courant--Friedrichs--Lewy condition, the time step is chose to be

@f[
\Delta t= \frac{\Delta x}{\lambda}  \quad\text{with}\quad
\lambda=\begin{cases}
1.2\,v_F&\text{if }S<0.36\,v_F\\
1.97\,S+v_F/2&\text{otherwise }
\end{cases}\label{eq:cfl}
@f]

such expression for the information speed @f$\lambda@f$ was empirically obtained for a wide variety of space discretizations and parameters.

@section twodrich Two-dimensional Richtmyer method


Regarding the two-dimensional case, one needs to consider the flux function along @f$y@f$
@f[
\frac{\partial}{\partial t}\vec{u}+ \frac{\partial}{\partial x}\vec{F}(\vec{u}) +\frac{\partial}{\partial y}\vec{G}(\vec{u})=0   
@f]
where now the state vector is @f$\vec{u}=[n,p_x,p_y]^\mathsf{T}@f$ and the fluxes functions are given by
@f[
\vec{F}(\vec{u})=\begin{bmatrix}      
p_xn^{-1/2}\\ \frac{p_x^2}{n^{3/2}}+\frac{v_F^2}{3}n^{3/2} + \frac{S^2}{2}n^2\\
p_xp_yn^{-3/2}
\end{bmatrix}
@f]
and
@f[
\vec{G}(\vec{u})=\begin{bmatrix}      
p_yn^{-1/2}\\
p_xp_yn^{-3/2} \\
\frac{p_y^2}{n^{3/2}}+\frac{v_F^2}{3}n^{3/2} + \frac{S^2}{2}n^2
\end{bmatrix}.
@f]

The generalisation of leads to a more complex method having, for the first step, the form
@f[
\vec{u}^{k+\sfrac{1}{2}}_{\subalign{i&+\sfrac{1}{2}\\j&+\sfrac{1}{2}}}=\frac{1}{4}\left(\vec{u}^{k}_{\subalign{i\\j}}+\vec{u}^{k}_{\subalign{&i+1\\&j}}+\vec{u}^{k}_{\subalign{&i\\&j+1}}+\vec{u}^{k}_{\subalign{i&+1\\j&+1}}\right)-\\
-\frac{\Delta t}{2\Delta x}\left(\vec{F}^{k}_{\subalign{i&+1\\j&+\sfrac{1}{2}}}-\vec{F}^{k}_{\subalign{&i\\&j+\sfrac{1}{2}}} \right)
-\frac{\Delta t}{2\Delta y}\left( \vec{G}^{k}_{\subalign{i&+\sfrac{1}{2}\\j&+1}}-\vec{G}^{k}_{\subalign{&i+\sfrac{1}{2}\\&j}} \right)
\label{eq:richt2D_1}
@f]
and then, the corrector is obtained by
@f[
\vec{u}^{k+1}_{\substack{i\\j}}=\vec{u}^{k}_{\subalign{&i\\&j}}-\frac{\Delta t}{\Delta x}\left(\vec{F}^{k+\sfrac{1}{2}}_{\subalign{&i+\sfrac{1}{2}\\&j}}-\vec{F}^{k+\sfrac{1}{2}}_{\subalign{&i-\sfrac{1}{2}\\&j}}\right)-\frac{\Delta t}{\Delta y}\left(\vec{G}^{k+\sfrac{1}{2}}_{\subalign{&i\\&j+\sfrac{1}{2}}}-\vec{G}^{k+\sfrac{1}{2}}_{\subalign{&i\\&j-\sfrac{1}{2}}}\right)\label{eq:richt2D_2}
@f]
In the matter of this scheme stability, for equal spacing @f$\Delta y=\Delta x@f$, the previous condition, defined for the one-dimensional case, has hitherto been sufficient to guarantee stability.


@section twodftcs Two-dimensional Forward Time Centred Space method

Turning now our attention to the diffusive part of which can be written as
@f[
\frac{\partial \vec{u}}{\partial t}=\nabla^2\vec{D}(\vec{u})+\vec{b}(\vec{u})
@f]
with the viscous term
@f[
\vec{D}(\vec{u})=\begin{bmatrix}      
0\\
(\nu_s\,p_x-\nu_o\,p_y)n^{-3/2}\\
(\nu_s\,p_y+\nu_o\,p_x)n^{-3/2}
\end{bmatrix}
@f]
and the magnetic source term
@f[
\vec{b}(\vec{u})=\begin{bmatrix}      
0\\
-\omega_c\, p_yn^{-1/2}\\
\hphantom{-}\omega_c\, p_xn^{-1/2}
\end{bmatrix},
@f]
applying a forward time and centred space stencil yields the scheme:
@f[
\vec{u}^{k+1}_{\substack{i\\j}}=\vec{u}^{k}_{\subalign{&i\\&j}}+\frac{\Delta t}{\Delta x^2}\left(\vec{D}^{k}_{\subalign{&i+1\\&j}} -2\vec{D}^{k}_{\subalign{&i\\&j}}+\vec{D}^{k}_{\subalign{&i-1\\&j}} \right)+\\+\frac{\Delta t}{\Delta y^2}\left(\vec{D}^{k}_{\subalign{&i\\&j+1}} -2\vec{D}^{k}_{\subalign{&i\\&j}}+\vec{D}^{k}_{\subalign{&i\\&j-1}} \right) +\Delta t\vec{b}^{k}_{\subalign{&i\\&j}}
@f]

It is important to notice that this method, being such a direct calculation has a strong constraint on the stability. For a purely diffusive system as @f$\partial_t u =\nu_s\,\partial^2_{xx} u@f$ it is well established that stability requires
@f[
2\nu_s \Delta t \leq \frac{\Delta x^2\Delta y^2}{\Delta x^2+ \Delta y^2},
@f]
which is somewhat narrower than the CFL condition for the hyperbolic part and imposes a maximum value of the @f$\nu_s@f$ parameter. Notwithstanding, given that the typical values for shear viscosity on graphene are so low this restriction is of little consequence. Yet, for the dispersive case of odd viscosity no general criterion is known.

@section twoddff Two-dimensional Du Fort--Frankel method

@f[
\vec{u}^{k+1}_{\substack{i\\j}}=\vec{u}^{k-1}_{\subalign{&i\\&j}}+\frac{2\Delta t}{\Delta x^2}\left[\vec{D}^{k}_{\subalign{&i+1\\&j}} -\left(\vec{D}^{k+1}_{\subalign{&i\\&j}} + \vec{D}^{k-1}_{\subalign{&i\\&j}} \right)+\vec{D}^{k}_{\subalign{&i-1\\&j}} \right]+\\+\frac{2\Delta t}{\Delta y^2}\left[\vec{D}^{k}_{\subalign{&i\\&j+1}} -\left(\vec{D}^{k+1}_{\subalign{&i\\&j}} + \vec{D}^{k-1}_{\subalign{&i\\&j}} \right)+\vec{D}^{k}_{\subalign{&i\\&j-1}} \right]
@f]

for a linear diffusion function @f$\partial_t u=\eta \nabla^2u@f$ can be simplified to
@f[  
\vec{u}^{k+1}_{\substack{i\\j}}=\frac{1-2\sigma-2\varsigma}{1+2\sigma+2\varsigma}\vec{u}^{k-1}_{\subalign{&i\\&j}}+\frac{2\sigma}{1+2\sigma+2\varsigma}\left[\vec{u}^{k}_{\subalign{&i+1\\&j}} +\vec{u}^{k}_{\subalign{&i-1\\&j}} \right]+\\+\frac{2\varsigma}{1+2\sigma+2\varsigma}\left[\vec{u}^{k}_{\subalign{&i\\&j+1}} +\vec{u}^{k}_{\subalign{&i\\&j-1}} \right]
@f]
with @f$\varsigma=\eta\Delta t / \Delta y^2@f$ and @f$\sigma=\eta\Delta t / \Delta x^2@f$. Evidently, for the equal grid @f$\Delta x=\Delta y@f$ case:
@f[
\vec{u}^{k+1}_{\substack{i\\j}}=\frac{1-4\sigma}{1+4\sigma}\vec{u}^{k-1}_{\subalign{&i\\&j}}+\frac{2\sigma}{1+4\sigma}\left(\vec{u}^{k}_{\subalign{&i+1\\&j}} +\vec{u}^{k}_{\subalign{&i-1\\&j}}+\vec{u}^{k}_{\subalign{&i\\&j+1}} +\vec{u}^{k}_{\subalign{&i\\&j-1}} \right)
@f]