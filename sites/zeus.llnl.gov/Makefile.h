# FLASH makefile definitions for the 64 bit Intel icc/ifort compiler on Linux
# Created for zeus.llnl.gov by nathan and cal.
#
# Fedora Core 3
# Intel Linux Fortran 9.0
# hdf5-1.6.5 (icc90)
# mpich-1.2.7 (icc90)

#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

MPI_PATH   = /usr/local
HDF4_PATH  =
HDF5_PATH  = /usr/global/tools/hdf5/installs/chaos_3_x86_64_ib/hdf5-1.6.5/parallel

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

FCOMP   = ${MPI_PATH}/bin/mpipathf90
CCOMP   = ${MPI_PATH}/bin/mpipathcc
CPPCOMP = ${MPI_PATH}/bin/mpipathCC
LINK    = ${MPI_PATH}/bin/mpipathf90 
 
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

FFLAGS_OPT   = -c -r8 -i8 -m64 -march=opteron -O3 -g
FFLAGS_DEBUG = -c -g -r8 -i4 -check bounds -check format -check output_conversion -warn all -real_size 64
FFLAGS_TEST  = -c -r8 -i4 -O2 -real_size 64


CFLAGS_OPT   = -c -O3 -r8 -i8 -m64 -march=opteron -g
CFLAGS_DEBUG = -c -g -debug extended -D_LARGEFILE64_SOURCE
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

LFLAGS_OPT   = -r8 -i8 -Ur -g -o
LFLAGS_DEBUG = -r8 -i4 -Vaxlib -g -o
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

LIB_HDF4 = -L $(HDF4_PATH)/lib -lmfhdf -ldf -lz -ljpeg
LIB_HDF5 =  $(HDF5_PATH)/lib/libhdf5.a -lz 
LIB_PNG  = -lpng

LIB_MPI     =
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
