#### IFORT ####
ifc=ifort
iflags=-fast -O3 -implicitnone -r8 -qopenmp -parallel -heap-arrays 200  -ansi-alias -qopt-jump-tables='large' -xcore-avx2
iflagsdbg=-traceback -O0 -warn all -implicitnone -r8 -check bounds -g -fpe0 -debug extended #-pg
iflagsprof=-traceback -O3 -implicitnone -r8 -g -debug inline-debug-info -parallel-source-info=2 -qopenmp -xcore-avx2 #-pg

#### GFORTRAN ####
gfc=gfortran-7 
gflags=-O2 -ffree-form -fimplicit-none -fopenmp -ffree-line-length-none #-fmax-stack-var-size=4096
gflagsdbg=-fbacktrace -O0 -ffree-form -Wall -Wextra -Wimplicit-interface -Wunused-parameter -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -fimplicit-none -fcheck=all -g -ffpe-trap=invalid,zero,overflow,underflow

objpath=./obj
resultspath=./Results


all:
	make run

init:
	mkdir -p Results
	mkdir -p obj

lib:
	reset
	@$(ifc) $(iflags) -c mymathlib.f90       -module $(objpath) -o $(objpath)/mymathlib.o
	@printf '%s' 'Compiled 1/4...'
	@$(ifc) $(iflags) -c classdef.f90        -module $(objpath) -o $(objpath)/classdef.o
	@printf '\r%s' 'Compiled 2/4...'
	@$(ifc) $(iflags) -c library.f90         -module $(objpath) -o $(objpath)/library.o
	@printf '\r%s' 'Compiled 3/4...'
	@$(ifc) $(iflags) -c postproc.f90        -module $(objpath) -o $(objpath)/postproc.o
	@printf '\r%s\n' 'Compiled 4/4...'

lib_dbg:
	reset
	@$(ifc) $(iflagsdbg) -c mymathlib.f90    -module $(objpath) -o $(objpath)/mymathlib.o
	@printf '%s' 'Compiled 1/4...'
	@$(ifc) $(iflagsdbg) -c classdef.f90     -module $(objpath) -o $(objpath)/classdef.o
	@printf '\r%s' 'Compiled 2/4...'
	@$(ifc) $(iflagsdbg) -c library.f90      -module $(objpath) -o $(objpath)/library.o
	@printf '\r%s' 'Compiled 3/4...'
	@$(ifc) $(iflagsdbg) -c postproc.f90     -module $(objpath) -o $(objpath)/postproc.o
	@printf '\r%s\n' 'Compiled 4/4...'

lib_prof:
	reset
	@$(ifc) $(iflagsprof) -c mymathlib.f90   -module $(objpath) -o $(objpath)/mymathlib.o
	@printf '%s' 'Compiled 1/4...'
	@$(ifc) $(iflagsprof) -c classdef.f90    -module $(objpath) -o $(objpath)/classdef.o
	@printf '\r%s' 'Compiled 2/4...'
	@$(ifc) $(iflagsprof) -c library.f90     -module $(objpath) -o $(objpath)/library.o
	@printf '\r%s' 'Compiled 3/4...'
	@$(ifc) $(iflagsprof) -c postproc.f90    -module $(objpath) -o $(objpath)/postproc.o
	@printf '\r%s\n' 'Compiled 4/4...'

run:
	reset
	make fileclean
	make init
	make lib
	@$(ifc) -I$(objpath) $(iflags) main.f90 $(objpath)/*.o -o main.out
	@time -f "	run time: %e" ./main.out

run_dbg:
	reset
	make fileclean
	make init
	make lib_dbg
	@$(ifc) -I$(objpath) $(iflagsdbg) main.f90 $(objpath)/*.o -o main.out
	@time -f "	run time: %e" ./main.out

run_prof:
	reset
	make fileclean
	make init
	make lib_prof
	@$(ifc) -I$(objpath) $(iflagsprof) main.f90 $(objpath)/*.o -o main.out

gridgen_dbg:
	reset
	@$(ifc) -I$(objpath) $(iflagsdbg) gridgen.f90 $(objpath)/*.o -o gridgen.out
	@./gridgen.out
	
gridgen:
	reset
	@$(ifc) -I$(objpath) $(iflags) gridgen.f90 $(objpath)/*.o -o gridgen.out
	@./gridgen.out
	
trial:
	reset
	make lib
	$(ifc) -I$(objpath) $(iflags) trial.f90 $(objpath)/*.o -o trial.out
	@./trial.out

trial_dbg:
	reset
	make lib_dbg
	$(ifc) -I$(objpath) $(iflagsdbg) trial.f90 $(objpath)/*.o -o trial.out
	@./trial.out

# Gfortran part
glib_dbg:
	reset
	@$(gfc) $(gflagsdbg) -c mymathlib.f90  -J$(objpath) -o $(objpath)/mymathlib.o
	@$(gfc) $(gflagsdbg) -c classdef.f90   -J$(objpath) -o $(objpath)/classdef.o
	@$(gfc) $(gflagsdbg) -c library.f90    -J$(objpath) -o $(objpath)/library.o
	@$(gfc) $(gflagsdbg) -c postproc.f90   -J$(objpath) -o $(objpath)/postproc.o
	@$(gfc) -c -I$(objpath) $(gflags) main.f90 $(objpath)/*.o 

glib:
	reset
	@$(gfc) $(gflags) -c mymathlib.f90     -J$(objpath) -o $(objpath)/mymathlib.o
	@$(gfc) $(gflags) -c classdef.f90      -J$(objpath) -o $(objpath)/classdef.o
	@$(gfc) $(gflags) -c library.f90       -J$(objpath) -o $(objpath)/library.o
	@$(gfc) $(gflags) -c postproc.f90      -J$(objpath) -o $(objpath)/postproc.o

grun:
	reset
	make fileclean
	make init
	make glib
	@$(gfc) -I$(objpath) $(gflags) main.f90 $(objpath)/*.o -o main.out
	@time -f "	Run time: %E" ./main.out
	#./main.out

grun_dbg:
	reset
	make fileclean
	make init
	make glib_dbg
	@$(gfc) -I$(objpath) $(gflagsdbg) main.f90 $(objpath)/*.o -o main.out
	@time -f "	Run time: %E" ./main.out
	#./main.out

gtrial:
	reset
	make init
	make glib
	$(gfc) -I$(objpath) $(gflags) trial.f90 $(objpath)/*.o -o trial.out
	@./trial.out

clean:
	-rm $(objpath)/*.o $(objpath)/*.mod *.out
	-rm visitlog.py 

fileclean:
	-rm $(resultspath)/*.tec
	-rm $(resultspath)/*.curve
	-rm $(resultspath)/*.dat

