# FLASH makefile definitions for the Intel ifc compilers on Linux

#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

HDF4_PATH  = /usr/local/hdf4
HDF5_PATH  = /opt/i386/hdf5/hdf5-1.6.5-intel9.1

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

FCOMP   = mpif90
CCOMP   = mpicc
CPPCOMP = mpiCC
LINK    = mpif90

# pre-processor flag
PP     = -D

#----------------------------------------------------------------------------
# Compilation flags
#
#  Three sets of compilation/linking flags are defined: one for optimized
#  code, one for testing, and one for debugging.  The default is to use the 
#  _OPT version.  Specifying -debug to setup will pick the _DEBUG version,
#  these should enable bounds checking.  Specifying _TEST is used for 
#  flash_test, and is set for quick code generation, and (sometimes) 
#  profiling.  The Makefile generated by setup will assign the generic token 
#  (ex. FFLAGS) to the proper set of flags (ex. FFLAGS_OPT).
#----------------------------------------------------------------------------

FFLAGS_OPT   =  -c -r8 -i4 -O2 -tpp7 -xN -mp1 -ftz -ip -u -assume byterecl
FFLAGS_DEBUG =  -c -r8 -i4 -g -assume byterecl
FFLAGS_TEST  =  -c -r8 -i4  -assume byterecl

CFLAGS_OPT   = -c -O2 -tpp7 -xN -mp1 -ip
CFLAGS_DEBUG = -c -g
CFLAGS_TEST  = -c -O2

CFLAGS_HDF5  = -I $(HDF5_PATH)/include
CFLAGS_NCMPI =

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT, 
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT   = -Vaxlib -O2 -tpp7 -xN -ip -o
LFLAGS_DEBUG = -r8 -i4 -Vaxlib -o
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

LIB_HDF4  = -L$(HDF4_PATH)/lib -lmfhdf -ldf -lz -ljpeg
LIB_HDF5  = -L $(HDF5_PATH)/lib -lhdf5 -lz

LIB_OPT   =  
LIB_DEBUG = 
LIB_TEST  =

LIB_NCMPI =
LIB_MPE   =
LIB_MPI   =

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



