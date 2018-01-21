#### IFORT ####
ifc=ifort
iflags=-fast -O3 -qopenmp -implicitnone -r8
#iflags=-traceback -O0 -warn all -implicitnone -r8 -check bounds -g -fpe0 #-pg

#### GFORTRAN ####
gfc=gfortran-7 
gflags=-O2 -ffree-form -fimplicit-none -fopenmp
#gflags=-fbacktrace -O0 -ffree-form -Wall -Wextra -Wimplicit-interface -Wunused-parameter -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -fimplicit-none -fcheck=all -g -ffpe-trap=invalid,zero
objpath=./obj
resultspath=./Results


init:
	mkdir -p Results
	mkdir -p obj

all:
	make run

lib:
	reset
	@$(ifc) $(iflags) -c mymathlib.f90          -module $(objpath) -o $(objpath)/mymathlib.o
	@$(ifc) $(iflags) -c vf_classdef.f90        -module $(objpath) -o $(objpath)/vf_classdef.o
	@$(ifc) $(iflags) -c vr_classdef.f90        -module $(objpath) -o $(objpath)/vr_classdef.o
	@$(ifc) $(iflags) -c wingpanel_classdef.f90 -module $(objpath) -o $(objpath)/wingpanel_classdef.o
	@$(ifc) $(iflags) -c wakepanel_classdef.f90 -module $(objpath) -o $(objpath)/wakepanel_classdef.o
	@$(ifc) $(iflags) -c library.f90            -module $(objpath) -o $(objpath)/library.o
	@$(ifc) $(iflags) -c postproc.f90           -module $(objpath) -o $(objpath)/postproc.o

run:
	reset
	make fileclean
	make lib
	@$(ifc) -I$(objpath) $(iflags) main.f90 $(objpath)/*.o -o main.out
	@time -f "	Run time: %E" ./main.out
	#time ./main.out

trial:
	reset
	make lib
	$(ifc) -I$(objpath) $(iflags) trial.f90 $(objpath)/*.o -o trial.out
	@./trial.out


clean:
	-rm $(objpath)/*.o $(objpath)/*.mod *.out

fileclean:
	-rm $(resultspath)/*.tec
	-rm $(resultspath)/*.curve

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

gtrial:
	reset
	make glib
	$(gfc) -I$(objpath) $(gflags) trial.f90 $(objpath)/*.o -o trial.out
	@./trial.out
