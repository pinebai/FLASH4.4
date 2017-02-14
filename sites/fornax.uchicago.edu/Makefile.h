# FLASH makefile definitions for the 64 bit Intel icc/ifort compiler on Linux
#
# CentOS release 5.2 (Final)
# Intel(R) Fortran Compiler for applications running on Intel(R) 64, Version 10.1    Build 20080312 Package ID: l_fc_p_10.1.015
# Intel icc (ICC) 10.1 20080312
# hdf5-1.6.9 (icc90)
# mpich-1.2.7p1 (icc90)

#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

MPI_PATH   = /usr/local/mpich/1.2.7p1/intel
HDF4_PATH  =
HDF5_PATH  = /usr/local/hdf5/1.6.9/intel

ZLIB_PATH  =

PAPI_PATH  =
PAPI_FLAGS =

NCMPI_PATH = /usr/local/pnetcdf/1.0.3/intel
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

FFLAGS_OPT   = -c -r8 -i4 -O3 -real_size 64 #-unroll -align -prefetch -pad -ip
FFLAGS_DEBUG = -c -g -r8 -i4 -check bounds -check format -check output_conversion -warn all -check uninit -traceback -fp-stack-check -real_size 64
FFLAGS_TEST  = -c -r8 -i4 -O2 -real_size 64


CFLAGS_OPT   = -c -O3 -D_LARGEFILE64_SOURCE
CFLAGS_DEBUG = -c -g -debug extended -D_LARGEFILE64_SOURCE
CFLAGS_TEST  = -c -O2 -D_LARGEFILE64_SOURCE

CFLAGS_HDF5  = -I$(HDF5_PATH)/include
CFLAGS_NCMPI = -I$(NCMPI_PATH)/include
CFLAGS_MPI   = -I$(MPI_PATH)/include

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT, 
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT   = -r8 -i4 -Vaxlib -lsvml -Ur -o
LFLAGS_DEBUG = -r8 -i4 -Vaxlib -g -fp-stack-check -o
LFLAGS_TEST  = -r8 -i4 -Vaxlib -o

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

LIB_HDF4 = 
LIB_HDF5 = -L $(HDF5_PATH)/lib -lhdf5 -lz
LIB_PNG  = -lpng

LIB_MPI   = -L$(MPI_PATH)/lib -lfmpich -lmpich
LIB_NCMPI = -L$(NCMPI_PATH)/lib -lpnetcdf
LIB_MPE   = -L$(MPI_PATH)/lib -lmpe


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
