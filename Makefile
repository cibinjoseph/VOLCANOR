#### IFORT ####
ifc=ifort
iflags=-fast -O3 -implicitnone -r8 -qopenmp -parallel -heap-arrays 4096 -ansi-alias -qopt-jump-tables='large'
iflagsdbg=-traceback -O0 -warn all -implicitnone -r8 -check bounds -g -fpe0 #-pg

#### GFORTRAN ####
gfc=gfortran-7 
gflags=-O2 -ffree-form -fimplicit-none -fopenmp #-fmax-stack-var-size=4096
gflagsdbg=-fbacktrace -O0 -ffree-form -Wall -Wextra -Wimplicit-interface -Wunused-parameter -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -fimplicit-none -fcheck=all -g -ffpe-trap=invalid,zero

objpath=./obj
resultspath=./Results


all:
	make run

init:
	mkdir -p Results
	mkdir -p obj

lib:
	reset
	@$(ifc) $(iflags) -c mymathlib.f90          -module $(objpath) -o $(objpath)/mymathlib.o
	@$(ifc) $(iflags) -c vf_classdef.f90        -module $(objpath) -o $(objpath)/vf_classdef.o
	@$(ifc) $(iflags) -c vr_classdef.f90        -module $(objpath) -o $(objpath)/vr_classdef.o
	@$(ifc) $(iflags) -c wingpanel_classdef.f90 -module $(objpath) -o $(objpath)/wingpanel_classdef.o
	@$(ifc) $(iflags) -c wakepanel_classdef.f90 -module $(objpath) -o $(objpath)/wakepanel_classdef.o
	@$(ifc) $(iflags) -c library.f90            -module $(objpath) -o $(objpath)/library.o
	@$(ifc) $(iflags) -c postproc.f90           -module $(objpath) -o $(objpath)/postproc.o

lib_dbg:
	reset
	@$(ifc) $(iflagsdbg) -c mymathlib.f90          -module $(objpath) -o $(objpath)/mymathlib.o
	@$(ifc) $(iflagsdbg) -c vf_classdef.f90        -module $(objpath) -o $(objpath)/vf_classdef.o
	@$(ifc) $(iflagsdbg) -c vr_classdef.f90        -module $(objpath) -o $(objpath)/vr_classdef.o
	@$(ifc) $(iflagsdbg) -c wingpanel_classdef.f90 -module $(objpath) -o $(objpath)/wingpanel_classdef.o
	@$(ifc) $(iflagsdbg) -c wakepanel_classdef.f90 -module $(objpath) -o $(objpath)/wakepanel_classdef.o
	@$(ifc) $(iflagsdbg) -c library.f90            -module $(objpath) -o $(objpath)/library.o
	@$(ifc) $(iflagsdbg) -c postproc.f90           -module $(objpath) -o $(objpath)/postproc.o

run:
	reset
	make fileclean
	make lib
	@$(ifc) -I$(objpath) $(iflags) main.f90 $(objpath)/*.o -o main.out
	@time -f "	run time: %e" ./main.out
	#time ./main.out

run_dbg:
	reset
	make fileclean
	make lib_dbg
	@$(ifc) -I$(objpath) $(iflagsdbg) main.f90 $(objpath)/*.o -o main.out
	@time -f "	run time: %e" ./main.out
	#time ./main.out

trial:
	reset
	make lib
	$(ifc) -I$(objpath) $(iflags) trial.f90 $(objpath)/*.o -o trial.out
	@./trial.out

# Gfortran part
glib_dbg:
	reset
	@$(gfc) $(gflagsdbg) -c mymathlib.f90          -J$(objpath) -o $(objpath)/mymathlib.o
	@$(gfc) $(gflagsdbg) -c vf_classdef.f90        -J$(objpath) -o $(objpath)/vf_classdef.o
	@$(gfc) $(gflagsdbg) -c vr_classdef.f90        -J$(objpath) -o $(objpath)/vr_classdef.o
	@$(gfc) $(gflagsdbg) -c wingpanel_classdef.f90 -J$(objpath) -o $(objpath)/wingpanel_classdef.o
	@$(gfc) $(gflagsdbg) -c wakepanel_classdef.f90 -J$(objpath) -o $(objpath)/wakepanel_classdef.o
	@$(gfc) $(gflagsdbg) -c library.f90            -J$(objpath) -o $(objpath)/library.o
	@$(gfc) $(gflagsdbg) -c postproc.f90           -J$(objpath) -o $(objpath)/postproc.o

glib:
	reset
	@$(gfc) $(gflags) -c mymathlib.f90          -J$(objpath) -o $(objpath)/mymathlib.o
	@$(gfc) $(gflags) -c vf_classdef.f90        -J$(objpath) -o $(objpath)/vf_classdef.o
	@$(gfc) $(gflags) -c vr_classdef.f90        -J$(objpath) -o $(objpath)/vr_classdef.o
	@$(gfc) $(gflags) -c wingpanel_classdef.f90 -J$(objpath) -o $(objpath)/wingpanel_classdef.o
	@$(gfc) $(gflags) -c wakepanel_classdef.f90 -J$(objpath) -o $(objpath)/wakepanel_classdef.o
	@$(gfc) $(gflags) -c library.f90            -J$(objpath) -o $(objpath)/library.o
	@$(gfc) $(gflags) -c postproc.f90           -J$(objpath) -o $(objpath)/postproc.o

grun:
	reset
	make fileclean
	make glib
	@$(gfc) -I$(objpath) $(gflags) main.f90 $(objpath)/*.o -o main.out
	@time -f "	Run time: %E" ./main.out
	#./main.out

grun_dbg:
	reset
	make fileclean
	make glib_dbg
	@$(gfc) -I$(objpath) $(gflagsdbg) main.f90 $(objpath)/*.o -o main.out
	@time -f "	Run time: %E" ./main.out
	#./main.out
gtrial:
	reset
	make glib
	$(gfc) -I$(objpath) $(gflags) trial.f90 $(objpath)/*.o -o trial.out
	@./trial.out

clean:
	-rm $(objpath)/*.o $(objpath)/*.mod *.out
	-rn visitlog.py 

fileclean:
	-rm $(resultspath)/*.tec
	-rm $(resultspath)/*.curve

