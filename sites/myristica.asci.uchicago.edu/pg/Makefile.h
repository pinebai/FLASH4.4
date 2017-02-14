# FLASH makefile definitions for ix86 Linux (Portland Group compiler)

# note, in order to get pgprof to work properly, it was necessary to 
# download a new version from ftp://ftp.pgroup.com/ that plays nicely
# with GNOME.

#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------
HDF4_PATH = /usr
HDF5_PATH = /usr

ZLIB_PATH  =
MPI_PATH = /usr/local/mpich-pg/

PAPI_PATH  =
PAPI_FLAGS =

NCMPI_PATH =
MPE_PATH   =

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#   DO NOT USE the MPICH wrappers around the compilers, use the native compiler
#   and specify the location of MPI (compiled for that compiler) in MPI_PATH
#   This way we can compile the code against dummy_MPI when that becomes available.
#
#----------------------------------------------------------------------------
FCOMP   = pgf90
CCOMP   = pgcc
CPPCOMP = pgCC
LINK    = pgf90
 
# pre-processor flag
PP      = -D 

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

FFLAGS_OPT = -c -fast -r8 -i4  
FFLAGS_DEBUG = -g -c -r8 -i4 
FFLAGS_TEST = -c -r8 -i4 -fast -Mprof=lines

FFLAGS_MPI  = -I $(MPI_PATH)/include

F90FLAGS =

CFLAGS_OPT = -O2  -c
CFLAGS_DEBUG = -g -c
CFLAGS_TEST = -c 

# if we are using HDF5, we need to specify the path to the include files
CFLAGS_HDF5 = -I $(HDF5_PATH)/include
CFLAGS_MPI  = -I $(MPI_PATH)/include
CFLAGS_NCMPI =

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT, 
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT   = -o
LFLAGS_DEBUG = -g -o
LFLAGS_TEST  = -Mprof=lines -o


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

LIB_HDF4  = -L $(HDF4_PATH)/lib -lmfhdf -ldf -ljpeg -lz
LIB_HDF5  = -L $(HDF5_PATH)/lib -lhdf5
LIB_MPI   = -L $(MPI_PATH)/lib -lmpich

LIB_PAPI  = 
LIB_MATH  = -ldfftw -ldrfftw

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



