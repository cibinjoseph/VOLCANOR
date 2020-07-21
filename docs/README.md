[![](https://img.shields.io/badge/status-under%20development-green.svg)]()  [![](https://img.shields.io/badge/Last%20Updated-June%202018-green.svg)]()  

![VOLCANOR](logo/VOLCANOR.png)

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

## Improvements (To be made)
### Stability Improvements
- Provision to provide precomputed initial solution

### Feature and Solution Improvements
- Account for pitching velocity
- Precompute trajectory
- Write out in binary format
- Prandtl-glauert, Karman-Tsien compressibility corrections etc.
- Check 25% of panel span inset of vortices create difference
- Include Trim algorithm

### Performance Improvements
- Implement recording to array before writing
- Implement free wake relaxation

## Contribution and Usage
This code is under active development and a lot of features--including a proper user-friendly interface--have yet to be added. Users are encouraged to go through the code if interested, and let the author know of issues and bugs, if found. However, be warned that most of the features are untested and unvalidated and the author offers no guarantee on the results obtained at this point in time.

## Authors
All code here was created by [Cibin Joseph](https://github.com/cibinjoseph) (cibinjoseph92@gmail.com).

## License
GNU General Public License v3.0
See [LICENSE](LICENSE) for full text.
