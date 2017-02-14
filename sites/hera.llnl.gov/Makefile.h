#       Makefile.h file for alc.llnl.gov
#
#	FLASH makefile definitions for Linux (Intel compiler)


#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------


HDF5path   = /usr/global/tools/hdf5/installs/chaos_3_x86_64/hdf5-1.6.5.intel/parallel
MPIpath    = /usr/lib/mpi/mpi_intel

GSL_PATH   = $(HOME)/gsl/linux-x86_64

PAPI_PATH  = /usr/local/tools/papi
PAPI_FLAGS = -c -I$(PAPI_PATH)/include

ZLIB_PATH  =

#NCMPI_PATH =
MPE_PATH   =

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#  We use the f90 compiler as the linker, so some C libraries may explicitly
#  need to be added into the link line.
#----------------------------------------------------------------------------

FCOMP      = mpiifort
CCOMP      = mpiicc
CPPCOMP    = mpiicc
LINK       = mpiifort

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

FFLAGS_OPT   =  -c -r8 -i4 -O3  -tpp7 -ip  -pad -unroll -fno-alias -safe_cray_ptr \
                -prefetch -D_LARGEFILE64_SOURCE -auto
FFLAGS_DEBUG =  -c -g -r8 -i4 -O0
FFLAGS_TEST  =  -c -r8 -i4 

FFLAGS_PAPI  = -I$(PAPI_PATH)/include
FFLAGS_MPI   = -I$(MPIpath)/include

F90FLAGS     =

CFLAGS       = -I$(HDF5path)/include -c -r8 -i4  -O3 -tpp7 -march=pentium4 -mcpu=pentium4 -ip -unroll -D_LARGEFILE64_SOURCE -I../lib/gsl/include
CFLAGS_HDF5  = -I$(HDF5path)/include
CFLAGS_MPI   = -I$(MPIpath)/include
CFLAGS_GSL   = -I$(GSL_PATH)/include

#CFLAGS_NCMPI =

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT, 
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS       =  -r8 -i4 -lm -Vaxlib -lsvml -o #-Wl,--allow-multiple-definition

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

#LIB_HDF4   = -L$(HDF4path) -lmfhdf -ldf -ljpeg -lz
LIB_HDF5    = -L$(HDF5path)/lib -lhdf5 -lz -L$(MPIpath)/lib -lpthread \
-Wl,-rpath,$(HDF5path)/lib -lelf

LIB_PAPI    = $(PAPI_PATH)/lib/libpapi.a $(PAPI_PATH)/lib/_fixunssfdi.o
LIB_PNG     = -lpng -lz

LIB_OPT     =
LIB_GSL     = -L$(GSL_PATH)/lib -lgsl -lgslcblas
LIB_DEBUG   = -Wl,-rpath,/usr/local/tv/dflt/linux-x86/lib -L/usr/local/tv/dflt/linux-x86/lib
LIB_TEST    =

#LIB_NCMPI   =
LIB_MPE     =

LIB_MPI     =
MPI_LIB     =

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



