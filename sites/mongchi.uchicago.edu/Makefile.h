# FLASH makefile definitions for the 64 bit Intel icc/ifort compiler on Linux
#
# Fedora Core 5
# Intel Linux Fortran 9.1.039
# Intel Linux C       9.1.042
# hdf5-1.6.5 (icc90)
# mpich-1.2.7 (icc90)

#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

#MPI_PATH   = /usr/local/mpich-1.2.7p1/intel
MPI_PATH   = /opt/pkgs/mpich2-1.3.1
HDF4_PATH  =
HDF5_PATH  = /opt/pkgs/hdf5-1.6.10

ZLIB_PATH  =

PAPI_PATH  =
PAPI_FLAGS =

FISHPAK_PATH = 

NCMPI_PATH = /usr/local/pnetcdf-1.0.1/intel2
MPE_PATH   =
HYPRE_PATH = $(HOME)/hypre

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#  We use the f90 compiler as the linker, so some C libraries may explicitly
#  need to be added into the link line.
#----------------------------------------------------------------------------

FCOMP   = ${MPI_PATH}/bin/mpif90 -fpic -i-dynamic -mcmodel=large -fpe0
CCOMP   = ${MPI_PATH}/bin/mpicc
CPPCOMP = ${MPI_PATH}/bin/mpiCC
LINK    = ${MPI_PATH}/bin/mpif90 -fpic -i-dynamic -mcmodel=large
 
# pre-processor flag

PP     = -D


#----------------------------------------------------------------------------
# GSL library for John ZuHone
GSL_PATH   = /usr/local/gsl-1.13
CFLAGS_GSL = -I${GSL_PATH}/include
LIB_GSL    = -L${GSL_PATH}/lib -lgsl -lgslcblas
#----------------------------------------------------------------------------



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

FFLAGS_OPT   = -c -r8 -i4 -O4 -real_size 64 #-check uninit -check bounds # -unroll -align -prefetch -pad -ip

#FFLAGS_DEBUG = -c -g -r8 -i4 -check bounds -check format -check output_conversion \
#-warn all -real_size 64 -check uninit -traceback -fp-stack-check -mcmodel=large

FFLAGS_DEBUG = -c -g -r8 -i4 -check bounds -check format -check output_conversion \
-real_size 64 -check uninit -traceback -fp-stack-check -mcmodel=large


#FFLAGS_DEBUG = -g -c -fdefault-real-8 -fdefault-double-8 \
#-ffree-line-length-none -pedantic -Wall -Wextra -Wconversion -Wunderflow \
#-ffpe-trap=invalid,zero,overflow -fbounds-check -fopenmp

FFLAGS_TEST  = -c -r8 -i4 -O2 -real_size 64
FFLAGS_HYPRE = -I${HYPRE_PATH}/include
CFLAGS_HYPRE = -I${HYPRE_PATH}/include

CFLAGS_OPT   = -c -O4 -D_LARGEFILE64_SOURCE
CFLAGS_DEBUG = -c -g -debug extended -D_LARGEFILE64_SOURCE -mcmodel=large
CFLAGS_TEST  = -c -O2 -D_LARGEFILE64_SOURCE

CFLAGS_HDF5  = -I $(HDF5_PATH)/include
CFLAGS_NCMPI = -I $(NCMPI_PATH)/include
CFLAGS_MPI   = -I$(MPI_PATH)/include

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT, 
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

#LFLAGS_OPT   = -r8 -i4 -Vaxlib -lsvml -Ur -o
#LFLAGS_DEBUG = -r8 -i4 -Vaxlib -g -mcmodel=large -o
#LFLAGS_TEST  = -r8 -i4 -Vaxlib -o

LFLAGS_OPT   = -diag-disable 10120 -O3 -o
LFLAGS_DEBUG = -diag-disable 10120 -o
LFLAGS_TEST  = -diag-disable 10120 -O2 -o


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
LIB_HDF5 = -L$(HDF5_PATH)/lib -lhdf5 -lz
LIB_PNG  = -lpng

#LIB_MPI     = -L$(MPI_PATH)/lib -lfmpich -lmpich -lmpichf90 -lmpichf90nc
#LIB_NCMPI = -L$(NCMPI_PATH)/lib -lpnetcdf
#LIB_MPE   = -L$(MPI_PATH)/lib -lmpe

LIB_MPI   =
LIB_NCMPI = -L$(NCMPI_PATH)/lib -lpnetcdf
LIB_MPE   = -L$(MPI_PATH)/lib -lmpe
LIB_HYPRE = -L${HYPRE_PATH}/lib -lHYPRE


LIB_FISHPAK = 

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
