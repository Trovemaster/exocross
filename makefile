goal:   xsec.x

tarball:
	tar cf duo.tar makefile *.f90
        
checkin:
	ci -l Makefile *.f90


FPATH = 

PLAT = _0501
FOR  = ifort

FFLAGS =  -O3 -ip -openmp 

#FFLAGS = -O0 -fpe0  -fltconsistency -stand f03 -check all -warn all -traceback -fp-stack-check  # debugging options


LAPACK = -static

LDIR = /nfs/workspaces/exomol/syurchenko/programs/libcerf/libcerf-1.4/lib/.lib/
LIB     =  # w_of_z.o err_fcts.o erfcx.o im_w_of_x.o   #  $(LIBDIR)libcerf.so.1.0.4 -L $LDIR 

###############################################################################

OBJ =  accuracy.o  timer.o input.o spectrum.o # use_libcerf_mod.o

xsec.x:	$(OBJ) crosssections.o 
	$(FOR) -o j-xsec$(PLAT).x $(OBJ) $(FFLAGS) crosssections.o $(LIB) -static

crosssections.o:	crosssections.f90 $(OBJ) 
	$(FOR) -c crosssections.f90 $(FFLAGS)

spectrum.o:	spectrum.f90 accuracy.o input.o  use_libcerf_mod.o
	$(FOR) -c spectrum.f90 $(FFLAGS)

accuracy.o:  accuracy.f90
	$(FOR) -c accuracy.f90 $(FFLAGS)

timer.o:  timer.f90
	$(FOR) -c timer.f90 $(FFLAGS)

input.o:  input.f90
	$(FOR) -c input.f90 $(FFLAGS)

use_libcerf_mod.o:  use_libcerf_mod.f90
	$(FOR) -c use_libcerf_mod.f90 $(FFLAGS)


clean:
	rm $(OBJ) *.mod duo.o

