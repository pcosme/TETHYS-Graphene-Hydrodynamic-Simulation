# TETHYS - *Two-dimensional Emitter of THz, Hydrodynamic Simulation.*
## Version 2.2.1

[![CodeScore](https://www.code-inspector.com/project/1694/score/svg)](https://www.code-inspector.com/project/1694/score/svg)
[![CodeStatus](https://www.code-inspector.com/project/1694/status/svg)](https://www.code-inspector.com/project/1694/status/svg)
[![CodeFactor](https://www.codefactor.io/repository/github/pcosme/hydrodynamic-simulation/badge?s=00232ac4455fd1f2e77fbc58fe3823f751721249)](https://www.codefactor.io/repository/github/pcosme/hydrodynamic-simulation)

![.github/workflows/main.yml](https://github.com/pcosme/Hydrodynamic-Simulation/workflows/.github/workflows/main.yml/badge.svg)

## Richtmyer method implementation
Repository for the elaboration of the hydrodynamic model simulation. 
Implemented for in 1D+1 and 2D+1 for density and velocity fields.

## Simplified flowchart of the code

![Flowchart](/images/CodeFlowchart.png)

## Requirements 

* gcc compiler 

* HDF5 libraries 

  For the writing of the complete data.

* ~~FFTW libraries~~

  Temporarily not implemented, in future will be used for calculate Fourier Transform of the signal. 

* Gnuplot & Python

  Responsible for the plotting

## Class Hierarchy

![Classes](/images/class_diagram_fluid+bc.svg)

## Style guide

### Semantic Versioning

Standard form of numeric *major.minor.patch* starting with the initial commit 1.0.0. Small (but relevant) bugs are considered lower level patches and new features (such as updating physical model) are minor level. Major level versions should be saved for breaking updates (like 2D implementation or parallelization)

#### Version history
1.0.0 Initial commit

1.1.0 Updated nonlinear terms in velocity flux, explicit pressure term. (minor changes in variables).

1.2.0 Addition of collisional loss term.
  <br>&emsp;1.2.1 Mean free path added as command line parameter. 
  <br>&emsp;1.2.2 Update of frequency functions + minor changes in appearance. 
  <br>&emsp;1.2.3 Added warning for the case on non-propagating plasmons + beginning implementation of CI tasks.
  <br>&emsp;1.2.4 New CFL condition.
  <br>&emsp;1.2.5 Output in HDF5 format.
  
1.3.0 Transition to object oriented code. Addition of Time Series Analysis and Electronic properties extraction.
  <br>&emsp;1.3.1 Stored energy at the capacitor gate calculated. 
  <br>&emsp;1.3.2 Viscosity term for Reynolds >10 implemented. 
  <br>&emsp;1.3.3 New class hierarchy on 1D algorithms. 
  <br>&emsp;1.3.4 Boundary conditions implemented as a separate class


2.0.0 Two dimensional code implementation. 1D version maintained for fast/simpler simulations. 
  <br>&emsp;2.0.1 New class hierarchy on 2D algorithms. New organization of header files
  <br>&emsp;2.0.2 Boundary conditions implemented as a separate class
  <br>&emsp;2.0.3 Linear _for_ loops
  <br>&emsp;2.1.0 Magnetic Field inclusion with _Godunov Splitting_ 
  <br>&emsp;2.2.0 Shear viscosity with FTCS method. Variable Aspect ratio 
  <br>&emsp;2.2.1 Momentum relaxation in 2D simulations
 
### Internal syntax

| Type            | Style                                 | E.g.              |
| :-------------: |:-------------:                        | :-----            |
| *Macros*          | Prefix + _ + Upper-case                | MAT_PI            |
| *Functions*       | Camel case                            | InitialCondRand   |
| *Variables*       | Lower-case 3 letters code + suffix     | den_mid           |

