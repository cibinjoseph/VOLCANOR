#### IFORT ####
ifc=ifort
iflags=-fast -O3 -implicitnone -r8 -qopenmp -parallel -heap-arrays 200  -ansi-alias -qopt-jump-tables='large' -xcore-avx2
iflagsdbg=-traceback -O0 -warn all -implicitnone -r8 -check bounds -g -fpe0 -debug extended -heap-arrays 200 #-pg
iflagsprof=-traceback -O3 -implicitnone -r8 -g -debug inline-debug-info -parallel-source-info=2 -qopenmp -xcore-avx2 #-pg

#### GFORTRAN ####
gfc=gfortran-7 
gflags=-O2 -ffree-form -fimplicit-none -fopenmp -ffree-line-length-none #-fmax-stack-var-size=4096
gflagsdbg=-fbacktrace -O0 -ffree-form -Wall -Wextra -Wimplicit-interface -Wunused-parameter -Wcharacter-truncation -Wsurprising -Waliasing -fimplicit-none -fcheck=all -g -ffpe-trap=invalid,zero,overflow,underflow

objpath=./obj
resultspath=./Results


all:
	make run

init:
	mkdir -p Results
	mkdir -p obj

lib:
	reset
	@$(ifc) $(iflags) -c libMath.f90       -module $(objpath) -o $(objpath)/libMath.o
	@printf '%s' 'Compiled 1/5...'
	@$(ifc) $(iflags) -c libC81.f90       -module $(objpath) -o $(objpath)/libC81.o
	@printf '\r%s' 'Compiled 2/5...'
	@$(ifc) $(iflags) -c classdef.f90        -module $(objpath) -o $(objpath)/classdef.o
	@printf '\r%s' 'Compiled 3/5...'
	@$(ifc) $(iflags) -c libCommon.f90         -module $(objpath) -o $(objpath)/libCommon.o
	@printf '\r%s' 'Compiled 4/5...'
	@$(ifc) $(iflags) -c libPostprocess.f90        -module $(objpath) -o $(objpath)/libPostprocess.o
	@printf '\r%s\n' 'Compiled 5/5...'

lib_dbg:
	reset
	@$(ifc) $(iflagsdbg) -c libMath.f90    -module $(objpath) -o $(objpath)/libMath.o
	@printf '%s' 'Compiled 1/5...'
	@$(ifc) $(iflagsdbg) -c libC81.f90    -module $(objpath) -o $(objpath)/libC81.o
	@printf '\r%s' 'Compiled 2/5...'
	@$(ifc) $(iflagsdbg) -c classdef.f90     -module $(objpath) -o $(objpath)/classdef.o
	@printf '\r%s' 'Compiled 3/5...'
	@$(ifc) $(iflagsdbg) -c libCommon.f90      -module $(objpath) -o $(objpath)/libCommon.o
	@printf '\r%s' 'Compiled 4/5...'
	@$(ifc) $(iflagsdbg) -c libPostprocess.f90     -module $(objpath) -o $(objpath)/libPostprocess.o
	@printf '\r%s\n' 'Compiled 5/5...'

lib_prof:
	reset
	@$(ifc) $(iflagsprof) -c libMath.f90   -module $(objpath) -o $(objpath)/libMath.o
	@printf '%s' 'Compiled 1/5...'
	@$(ifc) $(iflagsprof) -c libC81.f90   -module $(objpath) -o $(objpath)/libC81.o
	@printf '\r%s' 'Compiled 2/5...'
	@$(ifc) $(iflagsprof) -c classdef.f90    -module $(objpath) -o $(objpath)/classdef.o
	@printf '\r%s' 'Compiled 3/5...'
	@$(ifc) $(iflagsprof) -c libCommon.f90     -module $(objpath) -o $(objpath)/libCommon.o
	@printf '\r%s' 'Compiled 4/5...'
	@$(ifc) $(iflagsprof) -c libPostprocess.f90    -module $(objpath) -o $(objpath)/libPostprocess.o
	@printf '\r%s\n' 'Compiled 5/5...'

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
	@$(gfc) $(gflagsdbg) -c libMath.f90  -J$(objpath) -o $(objpath)/libMath.o
	@$(gfc) $(gflagsdbg) -c libC81.f90  -J$(objpath) -o $(objpath)/libC81.o
	@$(gfc) $(gflagsdbg) -c classdef.f90   -J$(objpath) -o $(objpath)/classdef.o
	@$(gfc) $(gflagsdbg) -c libCommon.f90    -J$(objpath) -o $(objpath)/libCommon.o
	@$(gfc) $(gflagsdbg) -c libPostprocess.f90   -J$(objpath) -o $(objpath)/libPostprocess.o
	@$(gfc) -c -I$(objpath) $(gflags) main.f90 $(objpath)/*.o 

glib:
	reset
	@$(gfc) $(gflags) -c libMath.f90     -J$(objpath) -o $(objpath)/libMath.o
	@$(gfc) $(gflags) -c libC81.f90     -J$(objpath) -o $(objpath)/libC81.o
	@$(gfc) $(gflags) -c classdef.f90      -J$(objpath) -o $(objpath)/classdef.o
	@$(gfc) $(gflags) -c libCommon.f90       -J$(objpath) -o $(objpath)/libCommon.o
	@$(gfc) $(gflags) -c libPostprocess.f90      -J$(objpath) -o $(objpath)/libPostprocess.o

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
	-rm -f $(objpath)/*.o $(objpath)/*.mod *.out
	-rm -f visitlog.py 

fileclean:
	-rm -f $(resultspath)/*.plt
	-rm -f $(resultspath)/*.curve
	-rm -f $(resultspath)/*.dat
	-rm -f $(resultspath)/*.txt
	-rm -f status.txt

