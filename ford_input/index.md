title: Introduction

The vortex lattice method is a potential flow based method. While different flavours exist, the one implemented here is a lifting surface model with a free-wake formulation. To keep the solver generalized, the flow solution is obtained in the global frame. Details on this specific implementation of the vortex lattice method are available in the following references:

1. Joseph, C., and Mohan, R. [A Parallel, Object-Oriented Framework for Unsteady Free-Wake Analysis of Multi-Rotor/Wing Systems.](https://doi.org/10.1016/j.compfluid.2020.104788) Computers &amp; Fluids, Elsevier BV, 215, Jan, 2021, p. 104788.  
2. Katz, J., and Plotkin, A. [Low-Speed Aerodynamics](https://doi.org/10.1017/CBO9780511810329). Cambridge University Press, Feb 05, 2001.

### Assumptions and limitations
1. Flow solution is incompressible, inviscid and irrotational. Viscous effects in flow like dissipation and strain are included using empirical factors. Reynolds number and Mach number effects on airfoil characteristics have not been modelled directly in the solver.
2. The Vatistas' core model is used for vortex desingularization. All vortices are rectilinear.
3. Airfoil thickness effects are not accounted for as only the camberline surface is modelled.
4. Input linear velocities are in inertial frame and angular velocities are in body frame. This is subject to change in future versions.

### Code philosophy
The code is centered around the object-oriented philosophy and is built up using a hierarchical model as shown in the figure. Vortex filaments are used to contstruct wing and wake ring elements. These are further used to construct wing and wake 'sheets' which make up larger aircraft configurations. Ever geometry is assumed to be a rotor at first. For example, a fixed-wing is a non-rotating single-bladed rotor. This rather unorthodox terminology was adopted for the sake of generalization during the initial phases of development.

<img src="|media|/data_abstraction.png" alt="" width="550pt">

### Code organization
The Fortran source files are present in the `src/` folder.

1. `main.f90` drives the solver and general program flow.  
2. `classdef.f90` contains attribute and method definitions for the various abstract data types.  
3. `libCommon.f90` contains subroutines that work on higher-level objects and those that deal with general bookkeeping.  
4. `libPostProcess.f90` contains post-processing and data write-out subroutines.  
5. `libMath.f90` is an independent math library that contains operations on vectors and matrices.  
6. `libC81.f90` deals with handling C81 airfoil files and is an independent library.  
7. `gridgen.90` converts results from the Lagrangian to an Eulerian framework to better plot pressure and velocity fields on a traditional CFD-like domain.

Data input to the solver is through the files `config.nml` and `geom01.nml` which are in _namelist_ format. Namelist files are a Fortran specific ASCII input format that follow a variable name-value ordering. Use the provided script `newcase.sh` to generate templates of these files.

A few utility codes are provided in `tools/` to make plotting and parsing of results easier.

### Installation
The solver has the following dependencies:

* CMake
* OpenMP
* LAPACK
* BLAS
* GNU or Intel Fortran compiler
* Python (Optional; for postprocessing)
* Paraview (Optional; for postprocessing)

The source files utilize a CMake build system. Follow the steps below to generate the `volcanor` executable.

1. Create a `build/` directory in the topmost folder.
2. Change directory to `build/`.
3. Run `cmake` from this directory pointing to the directory one level below. This will generate a Makefile for your OS and available compiler.
4. Run `make volcanor` to compile all source files and generate the `volcanor` executable in the `bin/` folder.

On a Linux system with Intel Fortran installed, this is achieved using the following commands.
```
mkdir build
cd build/
FC=ifort cmake ..
make volcanor
```
