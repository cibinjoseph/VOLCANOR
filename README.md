# Documentation
[![](https://img.shields.io/badge/status-under%20development-green.svg)]()  [![](https://img.shields.io/badge/Last%20Updated-Mar%202018-green.svg)]()  

A Parallel, Object oriented implementation of the Unsteady Vortex Lattice method in Fortran 90+ for aerodynamic analysis of a single wing under generic 3D motion.

## Features
- Parallelized implementation using OpenMP
- Provisions for analysis of sinusoid pitching and plunging motions
- Rankine vortex model (includes ideal vortex model)
- Vortex dissipation due to turbulence and viscosity
- Slow-start to avoid large starting vortex
- Wake strain to prevent violation of Helmholtz's Law
- Predictor-Corrector based wake convection for improved stability and accuracy
- Visualization of load history, circulation, wake structure etc in VisIt

## Known Issues
1. Overprediction of induced drag
2. Pitch rotation of blade about LE, should be customizable
3. Unequal spacing of panels causes instability when CP falls inside viscous core region

## Improvements
### Stability Improvements
- Implement CB2D
- Correct large starting vortex and root upwash
- Interpolation of vortex core radii along span for wake

### Feature and Solution Improvements
- Provision to restart
- Write out in binary format
- Prandtl-glauert, Karman-Tsien compressibility corrections etc.
- Induced drag computation drastically overpredicted
- Check 25% of panel span inset of vortices create difference
- Convert switches.f90 to readable input file

### Performance Improvements
- Deallocate unused variables depending on FD schemes
- Make wake array a shared variable in OpenMP part
- Implement recording to array before writing
- Implement free wake relaxation

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

## Authors
All code here was created by [Cibin Joseph](https://github.com/cibinjoseph) (cibinjoseph92@gmail.com).

## License
GNU General Public License v3.0
See [LICENSE](LICENSE) for full text.
