![VOLCANOR](logo/VOLCANOR.png)

# Description
A Parallel, Object oriented implementation of the Unsteady Vortex Lattice method in 
Fortran 90+ for aerodynamic analysis of rotors and wings under generic 3D motion.

## Usage

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

**rotor** (rotor_class)  
- nb - No. of blades   
- ns - No. of spanwise panels   
- nc - No. of chordwise panels   
- blade - blade objects  
- Omega - Angular velocity of rotor  
- Omega_slow -  
- shaft_axis -  
- hub_coords -  
- CG_coords -  
- radius -   
- chord -  
- root_cut -  
- control_pitch -  theta0,thetaC,thetaS  
- theta_twist -  
- pivotLE - pivot location from LE [x/c]  
- flap_hinge   hinge location from centre [x/R]  
- v_body -  
- om_body -  
- v_wind -  
- om_wind -  
- psi -  
- pts -  phi,theta,psi about CG_coords  
- spanwise_core -  
- streamwise_core -  
- AIC - Influence coefficient matrix  
- AIC_inv - Inverse of Influence coefficient matrix  
- gamvec -  
- gamvec_prev -  
- RHS -  
- init_wake_vel -  
- psi_start -  


**blade** (blade_class)   
- wiP - array of wing panels  
- waP - array of wake panels  
- Pwake - Predicted wake  
- theta - pitch angle  
- psi - azimuth angle  
- pivotLE - pivot distance from leading edge along chord  
- vind_wake - velocity induced at wake corners  
- vind_wake1, vind_wake2, vind_wake3 - velocity induced at wake corners for other schemes  
- Pvind_wake, vind_wake_step - velocity induced at wake corners for other schemes  
 
**wing panel**       (wingpanel_class)   
- vr - vortex ring
- pc - panel coords
- CP - collocation point coords
- ncap - unit normal vector
- velCP - local velocity at CP
- velCPm - relative inertial velocity at CP (due to motion)
- dForce - panel Force vector in inertial frame
- vel_pitch - pitch velocity
- dLift, dDrag - magnitudes of panel lift and drag
- delP - Pressure difference across panel
- panel_area - Panel area for computing lift
- r_hinge - dist to point about which pitching occurs (LE of wing)
- alpha - local angle of attack

**wake panel**       (wakepanel_class)   
- vr - vortex ring

**vortex ring**      (vr_class)   
- vf - fortex filaments (x4)
- gam - circulation

**vortex filament**  (vf_class)   
- fc - filament coordinates (xyz,1:2)  
- l0 - original length  
- lc - current length  
- r_vc0 - initial vortex core radius  
- r_vc  - current vortex core radius  
- age   - vortex age (in s)  

### Files

### Nomenclature

