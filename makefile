goal:   xsec.x

tarball:
	tar cf duo.tar makefile *.f90
        
checkin:
	ci -l Makefile *.f90


FPATH = 

PLAT = _1206_C
FOR  = ifort

FFLAGS = -C 
# -O3 -ip -openmp -static

LAPACK = # -mkl

LIB     =  

###############################################################################

OBJ =  accuracy.o  timer.o input.o spectrum.o 

xsec.x:	$(OBJ) crosssections.o 
	$(FOR) -o j-xsec$(PLAT).x $(OBJ) $(FFLAGS) crosssections.o $(LIB)

crosssections.o:	crosssections.f90 $(OBJ) 
	$(FOR) -c crosssections.f90 $(FFLAGS)

spectrum.o:	spectrum.f90 accuracy.o input.o
	$(FOR) -c spectrum.f90 $(FFLAGS)

accuracy.o:  accuracy.f90
	$(FOR) -c accuracy.f90 $(FFLAGS)

timer.o:  timer.f90
	$(FOR) -c timer.f90 $(FFLAGS)

input.o:  input.f90
	$(FOR) -c input.f90 $(FFLAGS)

clean:
	rm $(OBJ) *.mod duo.o

