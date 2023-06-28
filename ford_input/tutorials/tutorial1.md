title: Tutorial 1: Suddenly started rectangular wing

Let's simulate a rectangular wing, suddenly started from rest into a constant freestream velocity. 

### Case directory setup
Use the `newcase.sh` script to generate a new 'case' directory with the name `simplewing.case` as described in [Case directory setup](index.html).

### Case setup
The `config.nml` file contains global parameters that control the whole simulation like density, number of timesteps and timestep length. 

Geometry inside each simulation is defined using the `geomXX.nml` files where `XX` denotes a serial number.
