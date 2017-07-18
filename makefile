goal:   xsec.x

tarball:
	tar cf duo.tar makefile *.f90
        
checkin:
	ci -l Makefile *.f90


FPATH = 

PLAT = _0501
FOR  = ifort

CC   = icc
#CFLAGS = -O0 -g -std=c++0x

#FFLAGS =  -O0 -g 


CFLAGS = -O3 -xHost -std=c++0x -traceback -fp-model fast=2
 
FFLAGS =  -O3 -xHost -openmp -traceback

#FFLAGS = -O0 -g -fpe0  -fltconsistency -stand f03 -check all -warn all -traceback -fp-stack-check  # debugging options
#CFLAGS = -O0 -g -fpe0  -fltconsistency -traceback -fp-stack-check -std=c++0x

LAPACK = -static

#LDIR = /nfs/workspaces/exomol/syurchenko/programs/libcerf/libcerf-1.4/lib/.lib/
LIB     = -lstdc++ # w_of_z.o err_fcts.o erfcx.o im_w_of_x.o   #  $(LIBDIR)libcerf.so.1.0.4 -L $LDIR 

###############################################################################

OBJ =  accuracy.o  timer.o input.o spectrum.o VoigtKampff.o VoigtFortran.o profiles.o

xsec.x:	$(OBJ) crosssections.o 
	$(FOR) -o j-xsec$(PLAT).x $(OBJ) $(FFLAGS) crosssections.o $(LIB) -static

crosssections.o:	crosssections.f90 $(OBJ) 
	$(FOR) -c crosssections.f90 $(FFLAGS)

spectrum.o:	spectrum.f90 accuracy.o input.o  VoigtFortran.o
	$(FOR) -c spectrum.f90 $(FFLAGS)

accuracy.o:  accuracy.f90
	$(FOR) -c accuracy.f90 $(FFLAGS)

timer.o:  timer.f90
	$(FOR) -c timer.f90 $(FFLAGS)

input.o:  input.f90
	$(FOR) -c input.f90 $(FFLAGS)

profiles.o: VoigtKampff/profiles.cpp
	$(CC) -c VoigtKampff/profiles.cpp $(CFLAGS)

VoigtKampff.o: VoigtKampff/VoigtKampff.cpp profiles.o
	$(CC) -c VoigtKampff/VoigtKampff.cpp $(CFLAGS)

VoigtFortran.o: VoigtKampff/VoigtFortran.cpp VoigtKampff.o
	$(CC) -c VoigtKampff/VoigtFortran.cpp $(CFLAGS)
	

clean:
	rm $(OBJ) *.mod duo.o
