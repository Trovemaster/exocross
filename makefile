goal:   xsec.x

tarball:
	tar cf duo.tar makefile *.f90
        
checkin:
	ci -l Makefile *.f90


FPATH = 

PLAT = _2310_gf
FOR  = ifort

#FOR = gfortran
#FLAGS = -O2 -std=f2008 -fsyntax-only


FFLAGS =  -O3 -openmp -traceback  -ip 

#FFLAGS = -O0 -g -fpe0  -fltconsistency -stand f03 -check all -warn all -traceback -fp-stack-check  # debugging options

#FFLAGS = -O0 -fpe0  -fltconsistency -stand f03 -check all -warn all -traceback -fp-stack-check  # debugging options

#CFLAGS = -O0 -g -fpe0  -fltconsistency -traceback -fp-stack-check -std=c++0x

LAPACK = -static

FFLAGS2  =  #-vec-report3


#LDIR = /nfs/workspaces/exomol/syurchenko/programs/libcerf/libcerf-1.4/lib/.lib/

###############################################################################

OBJ =  accuracy.o  timer.o input.o spectrum.o VoigtKampff.o  phoenix.o

xsec.x:	$(OBJ) crosssections.o 
	$(FOR) -o j-xsec$(PLAT).x $(OBJ) $(FFLAGS) crosssections.o $(LIB) -static

crosssections.o:	crosssections.f90 $(OBJ) 
	$(FOR) -c crosssections.f90 $(FFLAGS)

spectrum.o:	spectrum.f90 accuracy.o input.o  VoigtKampff.o phoenix.o
	$(FOR) -c spectrum.f90 $(FFLAGS)

accuracy.o:  accuracy.f90
	$(FOR) -c accuracy.f90 $(FFLAGS)

timer.o:  timer.f90
	$(FOR) -c timer.f90 $(FFLAGS)

input.o:  input.f90
	$(FOR) -c input.f90 $(FFLAGS)

phoenix.o:  phoenix.f90 accuracy.o
	$(FOR) -c phoenix.f90 $(FFLAGS)

VoigtKampff.o:  VoigtKampff.f90 accuracy.o
	$(FOR) -c VoigtKampff.f90 $(FFLAGS)


clean:
	rm $(OBJ) *.mod duo.o
