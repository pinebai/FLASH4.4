HDF5_PATH = /usr/local/hdf5/1.6.9/intel
MPI_PATH = /usr/local/mpich/1.2.7p1/intel

CCOMP = mpicc
CFLAGS = -c -I$(HDF5_PATH)/include -g -I$(MPI_PATH)/include
LFLAGS =  -o
LIBS = -L$(HDF5_PATH)/lib -lhdf5 -lz -lmpich -lmpe

EXENAME = trajectorySort

OBJS = getOptions.o readFlashParticles.o sortParticles.o trajectory.o writeTrajectoryFile.o

.SUFFIXES: .c .o

.c.o : 
	$(CCOMP) $(CFLAGS) $(CDEFINES) $*.c

default: $(OBJS)
	$(CCOMP) $(LFLAGS) $(EXENAME) $(OBJS) $(LIBS) 

clean: 
	rm $(OBJS) $(EXENAME)
