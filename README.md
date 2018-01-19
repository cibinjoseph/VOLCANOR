# Documentation
Author : Cibin Joseph  
Last Updated : January 2018

## TO DO
- Validate pitching or plunging case (lift and drag)
- Induced drag computation drastically overpredicted
- Verify rotating wing results with BEMT
- Check whether file read write consumes large time 
- Parallelize all double do loops
- Implement recording to array before writing
- Check 25% of panel span inset of vortices create difference
- Check slow start of [3] extended tanh function
- Implement free wake relaxation
- Add rotor as xvec and yvec rotated about centre

## Rendering Readme.md
To render this document on a browser use command 
```sh
$ grip README.md
```
Print to pdf using usual Ctrl+P from browser.

## Code Details and Nomenclature
### Algorithm
*Note: All subroutines and functions calculate induced velocity by a unit vortex. Actual gamma has to be mutiplied.*
1. Wing coordinates in body frame
2. Rotate to inertial frame
3. TE vortex position
4. Influence coeff matrix
5. Initial Solution
   * Compute normal velocity at coll. points
   * Add kinematic velocity normal component at coll. point
   * Solve for gamma
6. LOOP START
7. Add position of shed vortex in wake array
8. Update coords of wing, vr, cp and ncap using UdelT
9. Compute induced velocities at cp
10. Solve for gamvec
11. Compute iduced velocity on wake vortices
12. Update vortex locations
13. LOOP END - GO TO 6.

### Features
1. Multiple vortex models 
2. Vortex dissipation
3. Slow-start to avoid large starting vortex
4. Wake strain

### Objects
#### vf_class  
- fc(3,2), lo, lc, r_vc
- vind(), calclength(), updlength()

#### vr_class
- vf(4), gam
- vind(), assignP, shiftdP, rot()

#### panel_class
- vr, pc(3,4), cp(3), ncap(3)
- assignP(), calcCP(), calcN(), rot()

### Parameters
Om = 600 rad/s => 62.832 rad/s => 3600 deg/s  
5 deg => 5/3600 s  
1 rev = 360 deg => 1/10 s  

