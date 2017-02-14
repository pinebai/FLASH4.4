# FLASH makefile definitions for ranger -- Sun constellation cluster at TACC
#                            for Portland Group compiler
#
# works with parallel hdf5
# optimizations not carefully tested, used defaults from generic PG compiler makefile
# 
# must add required modules before setup/build:
#   module remove gcc
#   module add phdf5
#   module add python
# (with recent defaults you have to remove the gcc module as shown here to get access to
# the phdf5 module.)
#
#   portland group is the default compiler/mpi stack for ranger
#   tested with  pgi/7.2-5  and  mvapich/1.0.1
#
# also 'module add phdf5' needs to appear in job submission script along with
# a 'module remove gcc' if necessary
#   (see accompanying example_batch_script)
#
#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

#HDF5_PATH =
#MPI_PATH  =

#SAMRAI_PATH =

#PAPI_PATH  =
#PAPI_FLAGS = -c -I$(PAPI_PATH)/include

#FISHPAK_PATH =
 

#ZLIB_PATH  =

#NCMPI_PATH =
#MPE_PATH   =

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#  We use the f90 compiler as the linker, so some C libraries may explicitly
#  need to be added into the link line.
#----------------------------------------------------------------------------

FCOMP   =  mpif90
CCOMP   =  mpicc
CPPCOMP =  mpicxx
LINK    =  mpif90


# Any additional machine specific libraries to be linked in
CONFIG_LIB =

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


#A note about Portland Group:
#the flags -fastsse and -pc 64 are used to avoid an 80-bit floating point 
#promotion that is a holdover from the x87 math co-processor instruction set

FFLAGS_OPT   = -c -r8 -i4 -O3
FFLAGS_DEBUG = -c -r8 -i4 -g -Ktrap=divz -Mchkfpstk -Mchkptr -Mchkstk -Mbounds
FFLAGS_TEST  = -c -r8 -i4

CFLAGS_OPT   = -c -O3
CFLAGS_DEBUG = -c -g
CFLAGS_TEST  = -c 

#FFLAGS_PAPI  = -I$(PAPI_PATH)/include
#FFLAGS_MPI   = -I$(MPI_PATH)/include

CFLAGS_HDF5  = -I${TACC_HDF5_INC} -DH5_USE_16_API
#CFLAGS_NCMPI = -I$(NCMPI_PATH)/include
#CFLAGS_MPI   = -I$(MPI_PATH)/include

#CFLAGS_SAMRAI = -I$(SAMRAI_PATH)/include

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT, 
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT   = -o
LFLAGS_DEBUG = -Wl,-noinhibit-exec -g -o
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


LIB_HDF5    = -L${TACC_HDF5_LIB} -lhdf5 -lz
              
LIB_MPI     =
LIB_PAPI    =
LIB_PNG     = -lpng -lz

LIB_OPT     = 
LIB_DEBUG   =
LIB_TEST    =

#LIB_SAMRAI  = -L$(SAMRAI_PATH)/lib -lSAMRAI

#LIB_FISHPAK = -L$(FISHPAK_PATH)/lib -lfishpak
#LIB_NCMPI   = -L$(NCMPI_PATH)/lib -lpnetcdf
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
AWK   = awk
CAT   = cat


