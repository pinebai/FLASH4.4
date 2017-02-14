# FLASH makefile definitions for the 64 bit gcc/gfortran compiler on Linux
#
# Red Hat Enterprise Linux 5.8
# Gnu Fortran 4.1.2
# Gnu C       4.1.2
# hdf5-1.8.8 (gcc412)
# Platform MPI (formerly HP MPI) 8.1.1

#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

# Load modules for compiler and HDF5 version; these must match, for example:
# module load libPHDF5/gnu/4.1.2/1.8.7
# Need the following for the Intel MKL, which contains Lapack and Blas
# module load intel/mkl/10.3.6.233

MPI_PATH   = /opt/platform_mpi
#HDF4_PATH  = /home/scarf027/hdf
HDF5_PATH  = `dirname ${HDF5_LIBDIR}`
HYPRE_PATH = ${HOME}/hypre_gnu
LAPACK_PATH = /apps/intel/2011.6.233/composer_xe_2011_sp1.6.233/mkl
LAPACK_LIBS = -lmkl_gf_lp64 -lmkl_sequential -lmkl_core
ZLIB_PATH  =

PAPI_PATH  =
PAPI_FLAGS =

FISHPAK_PATH = 

NCMPI_PATH = /usr/local/pnetcdf-icc
MPE_PATH   =

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#  We use the f90 compiler as the linker, so some C libraries may explicitly
#  need to be added into the link line.
#
#  Note that we use the MPI compilers to compile and link everything; this
#  is much easier than trying to mimic their behaviour through using the
#  serial compilers with extra flags.
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

FFLAGS_OPT   = -c -fdefault-real-8 -O3 -I$(LAPACK_PATH)/include -I$(HYPRE_PATH)/include # -unroll -align -prefetch -pad -ip
FFLAGS_DEBUG = -c -g -fdefault-real-8 -check bounds -check format -check output_conversion -warn all -check uninit -traceback -fp-stack-check -I$(LAPACK_PATH)/include -I$(HYPRE_PATH)/include
FFLAGS_TEST  = -c -fdefault-real-8 -O2 -I$(LAPACK_PATH)/include -I$(HYPRE_PATH)/include

CFLAGS_OPT   = -c -O3 -D_LARGEFILE64_SOURCE -I$(LAPACK_PATH)/include -I$(HYPRE_PATH)/include
CFLAGS_DEBUG = -c -g -debug extended -D_LARGEFILE64_SOURCE -I$(LAPACK_PATH)/include -I$(HYPRE_PATH)/include
CFLAGS_TEST  = -c -O2 -D_LARGEFILE64_SOURCE -I$(LAPACK_PATH)/include -I$(HYPRE_PATH)/include

CFLAGS_HDF5  = -I$(HDF5_PATH)/include -DH5_USE_16_API
CFLAGS_NCMPI = -I$(NCMPI_PATH)/include
CFLAGS_MPI   = -I$(MPI_PATH)/include

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT, 
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT   = -fdefault-real-8 -Ur -o
LFLAGS_DEBUG = -fdefault-real-8 -g -o
LFLAGS_TEST  = -fdefault-real-8 -o

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

#LIB_HDF4 = -L$(HDF4_PATH)/lib -lmfhdf -ldf -lz -ljpeg
LIB_HDF5 = -L$(HDF5_PATH)/lib -Wl,-rpath,$(HDF5_PATH)/lib -lhdf5 -lz
LIB_PNG  = -lpng

# We use mpif90 as the linker, so no need to link MPI libraries explicitly.
LIB_MPI   = 
LIB_NCMPI = -L$(NCMPI_PATH)/lib -lpnetcdf
LIB_MPE   =
LIB_HYPRE = -L$(HYPRE_PATH)/lib -Wl,-rpath,$(HYPRE_PATH)/lib -lHYPRE -L$(LAPACK_PATH)/lib/intel64 -Wl,-rpath,$(LAPACK_PATH)/lib/intel64 $(LAPACK_LIBS)

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
