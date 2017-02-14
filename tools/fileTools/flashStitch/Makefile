CCOMP = /usr/local/mpich-pg/bin/mpicc
CFLAGS = -c -O2 -g 
LFLAGS = -Bstatic -o
LIBS =

#Available defines: NO_HDF4, NO_HDF5, NEED_MPI
#Comment this out if your hdf5 doesn't need mpi. Sfocu should run as a single
#processor job in any case.

CDEFINES = 
HDF5_PATH = /usr/local/hdf5-pg
#NCMPI_PATH = /usr/local/pnetcdf-1.0.0
#MPI_PATH = /usr/local/mpich-pg

CFLAGS += -I$(HDF5_PATH)/include 

# HDF4 version
#LIBS = -lmfhdf -ldf -ljpeg -lz -lhdf5 -lm

LIBS += -L$(HDF5_PATH)/lib -lhdf5 -lz -lm
#LIBS += -L$(NCMPI_PATH)/lib -lpnetcdf
#LIBS += -L$(MPI_PATH)/lib -lmpich

EXENAME = flashStitcher

# add flash_reader_hdf4.o to the following to support HDF4
OBJS = flashStitcher.o main.o options.o

.SUFFIXES: .c .o
.c.o :
	$(CCOMP) $(CFLAGS) $(CDEFINES) $*.c

default: $(OBJS)
	$(CCOMP) $(LFLAGS) $(EXENAME) $(OBJS) $(LIBS) 

clean:
	rm $(OBJS) $(EXENAME)
