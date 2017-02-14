#
#  Makefile for forcing generator.
#
#  Written by
#
#          Christoph Federrath
#

FCOMP = gfortran -fdefault-real-8
# FCOMP = openmpif90 -fdefault-real-8 -O3
# FCOMP = ifort -r8 -i4 -O3
# FCOMP = xlf90 -qrealsize=8 -qintsize=4 -O3 -qstrict

forcing_generator : forcing_generator.o
	$(FCOMP) -o $@ forcing_generator.o

.SUFFIXES: .F90

.F90.o:
	$(FCOMP) -c $*.F90

clean :
	rm -f *.o *.mod *~ forcing_generator
