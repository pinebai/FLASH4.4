# FLASH makefile definitions for "lupin" ix86 Linux (Intel compilers)
#


#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

HDF5path   = /usr/local/hdf5-1.6.4-icc90
MPI_PATH    = /usr/local/mpich-1.2.7-icc-ifc90

SAMRAI_PATH = /usr/local/SAMRAI-v2.0.0

PAPI_PATH  = /usr/local/tools/papi
PAPI_FLAGS = -c -I$(PAPI_PATH)/include

FISHPAK_PATH = /usr/local/

NCMPI_PATH = /usr/local/pnetcdf-1.0.0-icc
MPE_PATH   =
ZLIB_PATH  =

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#  We use the f90 compiler as the linker, so some C libraries may explicitly
#  need to be added into the link line.
#----------------------------------------------------------------------------


FCOMP      = ifort
CCOMP      = icc
CPPCOMP    = icpc
LINK       = ifort -Bstatic 

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

FFLAGS_OPT   = -c -r8 -i4 -O3 -real_size 64
FFLAGS_DEBUG = -c -r8 -g -i4 -O0 -real_size 64
FFLAGS_TEST  = -c -r8 -i4 -real_size 64

CFLAGS = -c -O3 -D_LARGEFILE64_SOURCE

FFLAGS_PAPI  = -I$(PAPI_PATH)/include
FFLAGS_MPI   = -I$(MPI_PATH)/include

# specify the path to the include files
CFLAGS_HDF5  = -I$(HDF5path)/include
CFLAGS_NCMPI = -I$(NCMPI_PATH)/include
CFLAGS_MPI   = -I$(MPI_PATH)/include
CFLAGS_SAMRAI = -I$(SAMRAI_PATH)/include

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT, 
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT   = -r8 -i4 -Vaxlib -lsvml -Ur -o
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


LIB_HDF5    = -L$(HDF5path)/lib -lhdf5 -lz

LIB_MPI     = -L$(MPI_PATH)/lib -lmpich
LIB_PAPI    = $(PAPI_PATH)/lib/libpapi.a $(PAPI_PATH)/lib/_fixunssfdi.o
LIB_PNG     = -lpng -lz

LIB_OPT     = 
LIB_DEBUG   =
LIB_TEST    =

# libraries required to be linked in for specific machine
# once we understand these in more detail, we can move them elsewhere
# lupin needs -lgcc_eh library to work with Samrai
#CONFIG_LIB = -lrt -lpthread -lgcc_eh
CONFIG_LIB = -lrt -lpthread


SAMRAI_LIBLIST = -lSAMRAI1d_algs -lSAMRAI1d_appu -lSAMRAI1d_geom -lSAMRAI1d_hier -lSAMRAI1d_math_special -lSAMRAI1d_math_std -lSAMRAI1d_mblk -lSAMRAI1d_mesh -lSAMRAI1d_pdat_special -lSAMRAI1d_pdat_std -lSAMRAI1d_solv -lSAMRAI1d_xfer -lSAMRAI2d_algs -lSAMRAI2d_appu -lSAMRAI2d_geom -lSAMRAI2d_hier -lSAMRAI2d_math_special -lSAMRAI2d_math_std -lSAMRAI2d_mblk -lSAMRAI2d_mesh -lSAMRAI2d_pdat_special -lSAMRAI2d_pdat_std -lSAMRAI2d_solv -lSAMRAI2d_xfer -lSAMRAI3d_algs -lSAMRAI3d_appu -lSAMRAI3d_geom -lSAMRAI3d_hier -lSAMRAI3d_math_special -lSAMRAI3d_math_std -lSAMRAI3d_mblk -lSAMRAI3d_mesh -lSAMRAI3d_pdat_special -lSAMRAI3d_pdat_std -lSAMRAI3d_solv -lSAMRAI3d_xfer -lSAMRAI





LIB_SAMRAI  = -L$(SAMRAI_PATH)/lib $(SAMRAI_LIBLIST) -L/usr/lib/gcc/i386-redhat-linux/3.4.2 -lstdc++


LIB_FISHPAK = -L$(FISHPAK_PATH)/lib -lfishpak
LIB_NCMPI   = -L$(NCMPI_PATH)/lib -lpnetcdf
LIB_MPE     =

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



