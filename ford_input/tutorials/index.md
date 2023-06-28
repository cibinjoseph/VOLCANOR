title: Tutorials
ordered_subpage: tutorial1.md

@warning
Tutorial pages are under construction!

A few tutorials are provided for using the solver. Use the index on the left to access these tutorials. The basic case directory setup described below is the same for setting up all cases.

These tutorials were created for a Linux system.

### Case directory setup
Use the `newcase.sh` script to generate a new 'case' directory with an arbitrary name `mysimulation.case`. The `.case` extension is a simple convention to denote a single simulation case and any preferred name may be used.
```bash
./newcash.sh mysimulation.case
```

This bash script provides an easy way to set up a dedicated directory for a simulation by creating template input files and links to executables that are necessary to run a simulation.

```bash
tree mysimulation.case/

mysimulation.case/
├── bin -> /path/to/VOLCANOR/bin
├── config.nml
├── geom01.nml
├── gridconfig.nml
├── Makefile -> /path/to/VOLCANOR/tools/Makefile
├── Restart
├── Results
└── tools -> /path/to/VOLCANOR/tools
```
In the case directory, besides the links, you will find namelist (`.nml`) files which serve as inputs to VOLCANOR. Namelist files are a Fortran specifc ASCII input format which may also be read into Python ([f90nml](https://pypi.org/project/f90nml/)) and other languages using external packages. The `Results` directory will contain the solution outputs when they are generated on running the solver respectively. The `Restart` folder holds binary records of the simulation and may be used to restart the simulation if necessary and users may ignore this directory.

