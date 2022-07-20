![VOLCANOR](logo/VOLCANOR.png)

# Description
A Parallel, Object oriented implementation of the Unsteady Vortex Lattice method in 
Fortran 90+ for aerodynamic analysis of rotors and wings under generic 3D motion.

## Usage  
1. Use Cmake to build executable _volcanor_.
2. Create a newcase folder using _newcase.sh_.
3. Set up input parameters in _config.nml_ and _geom01.nml_.
4. Run using the make command.

## Code Details 
### Algorithm
Note: All subroutines and functions calculate induced velocity by a *unit* vortex. 
Actual gamma has to be mutiplied to find correct magnitude of induced velocity.
1. Wing coordinates in body frame
2. Transform to inertial frame
3. Prescribe TE vortex position
4. Compute influence coeff matrix
5. Initial Solution
  * Compute normal velocity at coll. points
  * Add normal component of kinematic velocity at coll. point
  * Solve for gamma
6. LOOP START
7. Add position of shed vortex in wake array
8. Update coords of wing, vr, cp and ncap using UdelT
9. Compute induced velocities at collocation point
10. Solve for gamvec
11. Compute induced velocity on wake vortices
12. Update vortex locations
13. LOOP END - GO TO 6.

### Object Hierarchy 

Objects | Nomenclature  
--- | --- 
Rotor  | rotor  
Blade  | blade     
Far wake panel        | Fwake  
Near wake panel       | Nwake    
Wing panel       | wingpanel    
Vortex ring      | vr    
Vortex filament  | vf    

### File formats

_.nml_ - Input namelist file  
_.dat_ - Binary file  
_.csv_ - CSV file (tab delimited) with first row as header  
_.json_ - JSON file  
_.log_ - Log file  
_.tec_ - Tecplot file

### Code Organization

**Input files**  
All input files have the extension _.in_    
_config.nml_  -  Main configuration file containing general solver inputs  
_geomXX.nml_  -  Configuration file for each geometry containing geometric parameters  
_gridconfig.nml_  -  Configuration file for generating grid-based solutions during postprocessing  

**Code files**  
_init_file.f90_  -  Variable initialisations  
_main.f90_  -  Main code that controls execution  
_classdef.f90_  -  Class definitions  
_libCommon.f90_  -  Subroutines common to all classes  
_libPostproc.f90_  -  Postprocessing subroutines  
_libMath.f90_  -  Math subroutines  
_gridgen.f90_  -  Generates grid-based data using gridconfig.in and filamentsXXXXX.dat  

**Output files**  
_wingNwakeXXXXX.plt_  -  Vortex collocation points of far wake   
_filamentsXXXXX.dat_  -  Filament properties in binary format  
_rXXforceDim.csv  -  Rotor dimensional force values at corresponding timesteps  
_rXXforceNonDim.csv  -  Rotor non-dimensional force values at corresponding timesteps  
_rXXbXXforceDistXXXXX.csv  -  Rotor sectional force values at corresponding timesteps  
