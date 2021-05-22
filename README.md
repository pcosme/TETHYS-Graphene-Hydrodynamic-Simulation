# TETHYS - *Two-dimensional Emitter of THz, Hydrodynamic Simulation.*
## Version 2.5.0

[![DOI](https://zenodo.org/badge/208904575.svg)](https://zenodo.org/badge/latestdoi/208904575)

![GitHub](https://img.shields.io/github/license/pcosme/TETHYS-Graphene-Hydrodynamic-Simulation)
![](https://img.shields.io/github/languages/top/pcosme/TETHYS-Graphene-Hydrodynamic-Simulation)

[![CodeScore](https://www.code-inspector.com/project/1694/score/svg)](https://www.code-inspector.com/project/1694/score/svg)
[![CodeStatus](https://www.code-inspector.com/project/1694/status/svg)](https://www.code-inspector.com/project/1694/status/svg)
[![CodeFactor](https://www.codefactor.io/repository/github/pcosme/tethys-graphene-hydrodynamic-simulation/badge?s=13e31a6a03d6b485a3f259a2e963d8584e2b0054)](https://www.codefactor.io/repository/github/pcosme/tethys-graphene-hydrodynamic-simulation)
![.github/workflows/main.yml](https://github.com/pcosme/TETHYS-Graphene-Hydrodynamic-Simulation/workflows/.github/workflows/main.yml/badge.svg)


## Documentation


The full documentation can be generated with Doxygen simply running 
```console
$ doxygen Doxyfile
```
from the root directory

<!--- ## Richtmyer method implementation --->
<!--- Repository for the elaboration of the hydrodynamic model simulation. --->
<!--- Implemented for in 1D+1 and 2D+1 for density and velocity fields. --->

## Simplified flowchart of the code

![Flowchart](./images/FlowchartTETHYS_2D.png)
![Flowchart](./images/FlowchartTETHYS_ELEC_2D.png)

## Requirements 

* gcc compiler 

* cmake (version 3.16)

  <small>To generate the makefiles</small>

* OpenMP<sup><small>&reg;</small></sup>

  <small>Responsible for CPU parallelization.</small>

* HDF5<sup><small>&reg;</small></sup> libraries (version 1.8.20 or higher) 

  <small>For the writing of the complete data.</small>
  
_<small>As an example in Fedora systems it should be enough to install the packages: cmake; libgomp; libomp; hdf5; hdf5-devel; hdf5-static </small>_
  
  

### Optional 

* Doxygen 

  <small>For generating the documentation</small>

* Gnuplot & Python

  <small>For the plotting</small>
  
## Usage   
### Compilation
For the compilation of the source code a makefile can ge automatically generated with cmake. It is generally a sensible idea to create a empty directory to hold the executable files

```console
$ mkdir build
$ cd build
```
And then, for the compilation itself. 
```console
$ cmake ..
$ make all
```

### Running a simulation

To run a simulation one can simply invoke
```console
$ ./TETHYS_2D 
```
and the program will prompt the user to provide the necessary parameters for the simulation. Alternatively, the user can pass directly such parameters as command arguments.
The necessary inputs are:
* Sound velocity S (vel_snd)
* Fermi velocity v<sub>F</sub> (vel_fer)
* Collision, of the electrons with impurities or defects, frequency (col)
* Shear viscosity (vis)
* Cyclotron frequency (cyc)
* (save_mode)

and optionally:  

* Aspect ratio of the simulation grid x:y (aspect_ratio)
```console
$ ./TETHYS_2D vel_snd vel_fer col vis cyc save_mode aspect_ratio
```


## Class Hierarchy

![Classes](./images/UML_Class_Diagram.png)

## Style guide

### Semantic Versioning

Standard form of numeric *major.minor.patch* starting with the initial commit 1.0.0. Small (but relevant) bugs are considered lower level patches and new features (such as updating physical model) are minor level. Major level versions should be saved for breaking updates (like 2D implementation or parallelization)

#### Version history
**1.0.0** Initial commit *[16 Sep. 2019]*

**1.1.0** Updated nonlinear terms in velocity flux, explicit pressure term. (minor changes in variables).

**1.2.0** Addition of collisional loss term.
  <br>&emsp;**1.2.1** Mean free path added as command line parameter. 
  <br>&emsp;**1.2.2** Update of frequency functions + minor changes in appearance. 
  <br>&emsp;**1.2.3** Added warning for the case on non-propagating plasmons + beginning implementation of CI tasks. *[17 Nov. 2019]*
  <br>&emsp;**1.2.4** New CFL condition.
  <br>&emsp;**1.2.5** Output in HDF5 format. *[8 Jan. 2020]*
  
**1.3.0** Transition to object oriented code. Addition of Time Series Analysis and Electronic properties extraction. *[1 May. 2020]*
  <br>&emsp;**1.3.1** Stored energy at the capacitor gate calculated.  
  <br>&emsp;**1.3.2** Viscosity term for Reynolds >10 implemented. *[26 May. 2020]*
  <br>&emsp;**1.3.3** New class hierarchy on 1D algorithms. *[26 Jun. 2020]*
  <br>&emsp;**1.3.4** Boundary conditions implemented as a separate class *[28 Jun. 2020]*

**2.0.0** Two dimensional code implementation. 1D version maintained for fast/simpler simulations. *[22 Jun. 2020]* 
  <br>&emsp;**2.0.1** New class hierarchy on 2D algorithms. New organization of header files *[26 Jun. 2020]*
  <br>&emsp;**2.0.2** Boundary conditions implemented as a separate class *[28 Jun. 2020]*
  <br>&emsp;**2.0.3** Linear _for_ loops *[29 Jun. 2020]*
  <br>&emsp;**2.1.0** Magnetic Field inclusion with _Godunov Splitting_ *[22 Jul. 2020]* 
  <br>&emsp;**2.2.0** Shear viscosity with FTCS method. Variable Aspect ratio  *[20 Sep. 2020]*
  <br>&emsp;**2.2.1** Momentum relaxation in 2D simulations *[24 Oct. 2020]*
  <br>&emsp;**2.3.0** (1,9) Weighted explicit method for viscous terms *[17 Nov. 2020]*
  <br>&emsp;**2.3.1** Parallelization with OpenMP  *[25 Nov. 2020]*
  <br>&emsp;**2.4.0** Simulation with odd viscosity  *[12 Feb. 2021]*
  <br>&emsp;**2.5.0** Energy transport by conduction and convection  *[22 May. 2021]*
### Internal syntax

| Type            | Style                                 | E.g.              |
| :-------------: |:-------------:                        | :-----            |
| *Macros*          | Prefix + _ + Upper-case                | MAT_PI            |
| *Functions*       | Camel case                            | InitialCondRand   |
| *Variables*       | Lower-case 3 letters code + suffix     | den_mid           |

