HDF5_PATH = /usr/local/hdf5-1.6.7/pg
MPI_PATH = /usr/local/mpich-1.2.7p1/pg

CCOMP = mpicc
CFLAGS = -c -I$(MPI_PATH)/include -I$(HDF5_PATH)/include -g 
LFLAGS =  -o
LIBS = -lmpich -lmpe -L$(HDF5_PATH)/lib -lhdf5 -lz

EXENAME = partSort

OBJS = writeFlashParticles.o readFlashParticles.o sortParticles.o writeFlashParticlesHDF5.o



.SUFFIXES: .c .o

.c.o : 
	$(CCOMP) $(CFLAGS) $(CDEFINES) $*.c

default: $(OBJS)
	$(CCOMP) $(LFLAGS) $(EXENAME) $(OBJS) $(LIBS) 

clean: 
	rm $(OBJS) $(EXENAME)
