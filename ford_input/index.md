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
The code is centered around the object-oriented philosophy and is built up using a hierarchical model. 

![Hierarchical data abstraction model](|media|/data_abstraction.png)

### Code organization

### Installation
