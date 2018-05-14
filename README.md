[![](https://img.shields.io/badge/status-under%20development-green.svg)]()  [![](https://img.shields.io/badge/Last%20Updated-May%202018-green.svg)]()  

![VOLCANOR](src_README/VOLCANOR.png)

# Documentation
A Parallel, Object oriented implementation of the Unsteady Vortex Lattice method in Fortran 90+ for aerodynamic analysis of rotors and wings under generic 3D motion.

## Features
- Parallelized implementation using OpenMP
- Rankine vortex model (includes ideal vortex model)
- Vortex dissipation due to turbulence and viscosity
- Slow-start to avoid large starting vortex
- Wake strain to prevent violation of Helmholtz's Law
- Predictor-Corrector based wake convection for improved stability and accuracy
- Visualization of load history, circulation, wake structure etc in VisIt

## Known Issues
1. Overprediction of induced drag
2. Unequal spacing of panels causes instability when CP falls inside viscous core region

## Improvements (To be made)
### Stability Improvements
- Interpolation of vortex core radii along span for wake

### Feature and Solution Improvements
- Slow start should not cause initial solution to go to zero for fixed wing.   
- Account for pitching velocity
- Precompute trajectory
- Provision to restart
- Write out in binary format
- Prandtl-glauert, Karman-Tsien compressibility corrections etc.
- Induced drag computation drastically overpredicted
- Check 25% of panel span inset of vortices create difference
- Auto evaluated unit tests for a few standard test cases

### Performance Improvements
- Verify wake array is used as a shared variable in OpenMP part
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

## Contribution and Usage
This code is under active development and a lot of features--including a proper user-friendly interface--have yet to be added. Users are encouraged to go through the code if interested, and let the author know of issues and bugs, if found. However, be warned that most of the features are untested and unvalidated and the author offers no guarantee on the results obtained at this point in time.

## Authors
All code here was created by [Cibin Joseph](https://github.com/cibinjoseph) (cibinjoseph92@gmail.com).

## License
GNU General Public License v3.0
See [LICENSE](LICENSE) for full text.
