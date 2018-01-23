# Documentation
[![](https://img.shields.io/badge/status-under%20development-green.svg)]()  

| | |
|---|---|
|Author        | Cibin Joseph         |  
|Last Updated  | Jan 2018             |  
|Pre-requisites| Fortran 90+, VisIt   |  

A Parallel, Object oriented implementation of the Unsteady Vortex Lattice method in Fortran 90+ for aerodynamic analysis of a single wing under generic 3D motion.

## Features
- Parallelized implementation using OpenMP
- Provisions for analysis of sinusoid pitching and plunging motions
- Choice of multiple vortex models 
- Vortex dissipation due to turbulence and viscosity
- Slow-start to avoid large starting vortex
- Wake strain to prevent violation of Helmholtz's Law
- Predictor-Corrector based wake convection for improved stability and accuracy
- Visualization of load history, wake structure etc in VisIt

## Known Issues
1. Overprediction of induced drag

## TO DO
- Add predictor-corrector approach
- Induced drag computation drastically overpredicted
- Verify rotating wing results with BEMT
- Check whether file read write consumes large time 
- Parallelize all double do loops
- Implement recording to array before writing
- Check 25% of panel span inset of vortices create difference
- Check slow start of [3] extended tanh function
- Implement free wake relaxation
- Add rotor as xvec and yvec rotated about centre

## Code Details 
### Algorithm
Note: All subroutines and functions calculate induced velocity by a *unit* vortex. Actual gamma has to be mutiplied to find correct magnitude of induced velocity.
1. Wing coordinates in body frame
2. Transform to inertial frame
3. Prescribe TE vortex position
4. Compute influence coeff matrix
5. Initial Solution
   * Compute normal velocity at coll. points
   * Add kinematic velocity normal component at coll. point
   * Solve for gamma
6. LOOP START
7. Add position of shed vortex in wake array
8. Update coords of wing, vr, cp and ncap using UdelT
9. Compute induced velocities at collocation point
10. Solve for gamvec
11. Compute induced velocity on wake vortices
12. Update vortex locations
13. LOOP END - GO TO 6.

### Parameters for test case
Om = 600 rad/s => 62.832 rad/s => 3600 deg/s  
5 deg => 5/3600 s  
1 rev = 360 deg => 1/10 s  

