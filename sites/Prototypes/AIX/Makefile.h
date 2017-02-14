

#----------------------------------------------------------------------------
# Set the HDF/HDF5 and PAPI library paths 
# -- these need to be updated for your system
# If PAPI doesn't exist on your system, comment them out
#----------------------------------------------------------------------------
HDF5_PATH  = /usr/local/hdf5/hdf5-1.4.2/parallel


PAPI_PATH  = /usr/local
PAPI_FLAGS = -c -I$(PAPI_PATH)/include -qsuffix=f=F90:cpp=F90 -qfree

ZLIB_PATH  =

NCMPI_PATH =
MPE_PATH   =

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#   Use the MPICH wrappers around the compilers -- these will automatically
#   load the proper libraries and include files. 
#----------------------------------------------------------------------------

FCOMP   = mpxlf90_r
CCOMP   = mpcc_r
CPPCOMP = mpCC_r
LINK    = mpxlf90_r

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

FFLAGS_OPT   = -O3 -qintsize=4 -qrealsize=8 -cpp -c \
               -qxlf90=autodealloc \
               -qsuffix=cpp=F -qtune=auto  

# qhot causes problems with 2.0
# -qhot -qcache=auto

FFLAGS_DEBUG = -g -qintsize=4 -qrealsize=8 -cpp -c 
FFLAGS_TEST  =

F90FLAGS     = -qsuffix=f=F90:cpp=F90 -qfree
f90FLAGS     = -qsuffix=f=f90:cpp=F90 -qfree


# if we are using HDF5, we need to specify the path to the include files
CFLAGS_HDF5  = -I $(HDF5_PATH)/include/ -DNOUNDERSCORE

CFLAGS_OPT   = -c -O3 -qcache=auto -qtune=auto -DIBM
CFLAGS_DEBUG = -g -c -DIBM
CFLAGS_TEST  =

CFLAGS_HDF5  =
CFLAGS_NCMPI =

MDEFS = -WF,

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT, 
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

#  Linker flags for optimization
LFLAGS_OPT   = -bmaxdata:0x80000000 -o

#  Linker flags for debugging
LFLAGS_DEBUG = -bmaxdata:0x80000000 -o

LFLAGS_TEST  = -bmaxdata:0x80000000 -o

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


LIB_HDF5  = -L $(HDF5_PATH)/lib -lhdf5 -L /usr/local/lib -lz 
LIB_PAPI  = -L$(PAPI_PATH)/lib -lpapi -L/usr/lpp/pmtoolkit/lib -lpmapi
LIB_MATH  = -lessl

LIB_OPT   = 
LIB_DEBUG =
LIB_TEST  =

LIB_NCMPI =
LIB_MPE   =
LIB_MPI  =

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

MV   = mv -f
AR   = ar -r
RM   = rm -f
CD   = cd
RL   = ranlib
ECHO = echo
