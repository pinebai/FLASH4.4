# FLASH makefile definitions for x86_64 Linux
#
# Fedora Core 10
# hdf5-1.8.3
# mpich2-1.0.8p1

#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

HDF4_PATH  =
HDF5_PATH  = /usr/local/hdf5
MPI_PATH   = /usr/local/mpich2

ZLIB_PATH  =

PAPI_PATH  =
PAPI_FLAGS =

NCMPI_PATH =
MPE_PATH   =

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#  We use the f90 compiler as the linker, so some C libraries may explicitly
#  need to be added into the link line.
#----------------------------------------------------------------------------

FCOMP   = ${MPI_PATH}/bin/mpif90
CCOMP   = ${MPI_PATH}/bin/mpicc
CPPCOMP = ${MPI_PATH}/bin/mpiCC
LINK    = ${MPI_PATH}/bin/mpif90

# pre-processor flag

PP     = -D

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

FFLAGS_OPT = -c -O2 -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none
FFLAGS_DEBUG = -g -Wall -fbounds-check -c -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none
FFLAGS_TEST  = -c -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none

F90FLAGS =

CFLAGS_OPT = -O2 -c
CFLAGS_DEBUG = -g -Wall -fbounds-check -c
CFLAGS_TEST = -c 

# if we are using HDF5, we need to specify the path to the include files
CFLAGS_HDF5 = -I $(HDF5_PATH)/include -DH5_USE_16_API
CFLAGS_NCMPI =
CFLAGS_MPI = -I$(MPI_PATH)/include

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT, 
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT   = -o
LFLAGS_DEBUG = -g -o
LFLAGS_TEST  = -o

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

LIB_OPT   =
LIB_DEBUG =
LIB_TEST  =

LIB_HDF4 = -L $(HDF4_PATH)/lib -lmfhdf -ldf -lz -ljpeg
LIB_HDF5 = -L $(HDF5_PATH)/lib -lhdf5 -lz
LIB_PNG  = -lpng

LIB_MPI     = -L$(MPI_PATH)/lib -lfmpich -lmpich
LIB_NCMPI =
LIB_MPE   =

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

MV = mv -f
AR = ar -r
RM = rm -f
CD = cd
RL = ranlib
ECHO = echo
