# importing the flash configuration
include ../../Makefile.h

# Generated automatically from Makefile.in by configure.
MKLIB=ar -r

DMPILIB = ../object/libdmpi.a
DMPIINC = ../include


IFILES =  -I$(DMPIINC) 

FFLAGS = -O -r8 -i4 
F77FLAGS = -O -r8 -i4 

####### Rules
.SUFFIXES:

.SUFFIXES: .f .F .f90 .F90 .c .o

.f.o :
	$(FCOMP) $(F77FLAGS) $(IFILES)  -c $<
.F.o :
	$(FCOMP) $(F77FLAGS) $(IFILES)  -c $<
.f90.o :
	$(FCOMP) $(FFLAGS) $(IFILES)  -c $<
.F90.o :
	$(FCOMP) $(FFLAGS) $(IFILES)  -c $<
.c.o :
	$(CCOMP) $(CFLAGS) $(IFILES)  -c $<

OBJS =   mpi_dummyF.o mpi_dummyC.o 

$(DMPILIB): $(OBJS)
	( if test ! -d  ../object ; then mkdir ../object; fi)
	rm -f $(DMPILIB); $(MKLIB) $(DMPILIB) $(OBJS)

$(OBJS) : $(DMPIINC)/mpif.h $(DMPIINC)/mpi.h


clean:
	rm -f $(OBJS)

realclean: clean
	rm -rf *~ $(DMPILIB)


