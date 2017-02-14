#       Makefile.h file for alc.llnl.gov
#
#	FLASH makefile definitions for Linux (Intel compiler)


#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------


HDF5path   = /usr/local/hdf5-1.6.5/gcc
MPIpath    = /usr/local/mpich-1.2.7p1/lahey

PAPI_PATH  = 
PAPI_FLAGS = 

ZLIB_PATH  =

NCMPI_PATH = /usr/local/pnetcdf-1.0.0/gcc
MPE_PATH   =

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#  We use the f90 compiler as the linker, so some C libraries may explicitly
#  need to be added into the link line.
#----------------------------------------------------------------------------

FCOMP      = ${MPIpath}/bin/mpif90
CCOMP      = gcc
CPPCOMP    = gcc
LINK       = ${MPIpath}/bin/mpif90


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

FFLAGS_OPT   = -c --o2 --tp4 -CcdRR8  
FFLAGS_DEBUG = -c -g  --trace  --trap --chk[aes] -CcdRR8  
FFLAGS_TEST  = -c -CcdRR8 

FFLAGS_PAPI  = -I$(PAPI_PATH)/include

F90FLAGS     =

#CFLAGS       = -c -O3 -tpp7 -march=pentium4 -mcpu=pentium4 -ip -unroll \
#               -D_LARGEFILE64_SOURCE

CFLAGS_OPT =   -I$(MPIpath)/include -c -O2 -D_LARGEFILE64_SOURCE -g
CFLAGS_DEBUG = -I$(MPIpath)/include -c -g 
CFLAGS_TEST  = -I$(MPIpath)/include -c
CFLAGS_HDF5  = -I$(HDF5path)/include
CFLAGS_NCMPI = -I$(NCMPI_PATH)/include

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT, 
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS       =  -o  

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

LIB_MPI     =
LIB_HDF5    = -L$(HDF5path)/lib -lhdf5 -lz  

LIB_PAPI    = $(PAPI_PATH)/lib/libpapi.a $(PAPI_PATH)/lib/_fixunssfdi.o
LIB_PNG     = -lpng -lz

LIB_OPT     =
LIB_DEBUG   =
LIB_TEST    =

LIB_NCMPI   = -L$(NCMPI_PATH)/lib -lpnetcdf -L/usr/local/pgi-5.2.4/linux86/5.2/lib -lpgc

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



