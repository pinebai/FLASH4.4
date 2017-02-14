#       Makefile.h file for ellipse.uchicago.edu
#
#	FLASH makefile definitions for Linux (Intel compiler)


#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------


HDF5_PATH  = /soft/hdf5-1.6.7-icc-gm-r1
MPI_PATH   = /soft/mpich-gm-1.2.7p1.2-intel-r1

PAPI_PATH  = 
PAPI_FLAGS = 

ZLIB_PATH  =

NCMPI_PATH = /soft/parallel-netcdf-1.0.1-icc-gm-r1
MPE_PATH   =

# for John Zuhones DPMTA library, only on ellipse
#DPMTA_PATH = /scratch7/zuhone/dpmta

GM_PATH    = /soft/gm-2.1.26_Linux-r1
PVFS1_PATH = /soft/pvfs-1.6.3-r1
PVFS2_PATH = /soft/pvfs-2.6.2-r1


#----------------------------------------------------------------------------
# Compiler and linker commands
#
#  We use the f90 compiler as the linker, so some C libraries may explicitly
#  need to be added into the link line.
#----------------------------------------------------------------------------

FCOMP      = $(MPI_PATH)/bin/mpif90
CCOMP      = $(MPI_PATH)/bin/mpicc
CPPCOMP    = $(MPI_PATH)/bin/mpicc
LINK       = ifort -static
LINK       = ${MPI_PATH}/bin/mpif90 -static

# pre-processor flag

PP         = -D

#-----------------------------------------------------------------------------
# Compilation flags
#
#  Three sets of compilation/linking flags are defined: one for optimized code
#  code ("-opt"), one for debugging ("-debug"), and one for testing ("-test").
#  Passing these flags to the setup script will cause the value associated with
#  the corresponding keys (i.e. those ending in "_OPT", "_DEBUG", or "_TEST") to
#  be incorporated into the final Makefile. For example, passing "-opt" to the
#  setup script will cause the flags following "FFLAGS_OPT" to be assigned to
#  "FFLAGS" in the final Makefile. If none of these flags are passed, the default
#  behavior will match that of the "-opt" flag.
#  In general, "-opt" is meant to optimize compilation and linking. "-debug"
#  should enable runtime bounds checking, debugger symbols, and other compiler-
#  specific debugging options. "-test" is useful for testing different
#  combinations of compiler flags particular to your individual system.
#----------------------------------------------------------------------------

FFLAGS_OPT   =  -c -r8 -i4 -O3  -real_size 64 -auto
FFLAGS_DEBUG =  -c -g -r8 -i4 -C -real_size 64 -auto -traceback
FFLAGS_TEST  =  -c -r8 -i4 -real_size 64

FFLAGS_PAPI  = -I$(PAPI_PATH)/include
FFLAGS_MPI   = -I$(MPI_PATH)/include
FFLAGS_DPMTA = -I$(DPMTA_PATH)/include

F90FLAGS     =

CFLAGS_OPT   = -c -O3 -D_LARGEFILE64_SOURCE
CFLAGS_DEBUG   = -c -g -C -traceback -debug all -D_LARGEFILE64_SOURCE

CFLAGS_HDF5  = -I$(HDF5_PATH)/include
CFLAGS_MPI   = -I$(MPI_PATH)/include
CFLAGS_NCMPI = -I$(NCMPI_PATH)/include

CFLAGS_DPMTA = -I$(DPMTA_PATH)/include

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT, 
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS       = -r8 -i4 -lm -Vaxlib -lsvml -o 

#----------------------------------------------------------------------------
# Library specific linking
#
#  If a FLASH module has a 'LIBRARY xxx' line in its Config file, we need to
#  create a macro in this Makefile.h for LIB_xxx, which will be added to the
#  link line when FLASH is built.  This allows us to switch between different
#  (incompatible) libraries.  We also create a _OPT, _DEBUG, and _TEST
#  library macro to add any performance-minded libraries (like fast math),
#  depending on how FLASH was setup.
#----------------------------------------------------------------------------

LIB_HDF5    = -L$(HDF5_PATH)/lib -lhdf5 -lz  -L$(MPI_PATH)/lib -lmpich -lpthread
              
LIB_MPI     = -L$(MPI_PATH)/lib -lmpich
LIB_PAPI    = $(PAPI_PATH)/lib/libpapi.a $(PAPI_PATH)/lib/_fixunssfdi.o
LIB_PNG     = -lpng -lz

LIB_OPT     = -L$(PVFS1_PATH)/lib -lpvfs -L$(PVFS2_PATH)/lib -lpvfs2 -L$(GM_PATH)/lib -lgm
LIB_DEBUG   = -L$(PVFS1_PATH)/lib -lpvfs -L$(PVFS2_PATH)/lib -lpvfs2 -L$(GM_PATH)/lib -lgm
LIB_TEST    =

LIB_NCMPI   = -L$(NCMPI_PATH)/lib -lpnetcdf 
LIB_MPE     =

LIB_DPMTA   = -L$(DPMTA_PATH)/lib -ldpmta
#----------------------------------------------------------------------------
# Additional machine-dependent object files
#
#  Add any machine specific files here -- they will be compiled and linked
#  when FLASH is built.
#----------------------------------------------------------------------------

MACHOBJ = 

#----------------------------------------------------------------------------
# Additional commands
#---------------------------------------------------------------------------- 

MV    = mv -f
AR    = ar -r
RM    = rm -f
CD    = cd
RL    = ranlib
ECHO  = echo



