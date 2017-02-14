#	FLASH makefile definitions for ASCI Red machine (Intel)
#
#	2/28/99
#
#----------------------------------------------------------------------------

HDF5_PATH  = /usr/community/hdf5/hdf5-1_4_0
HDF4_PATH  = /usr/community/hdf

ZLIB_PATH  = /usr/community/hdf5/ZLIB

PAPI_PATH  =
PAPI_FLAGS =

NCMPI_PATH =
MPE_PATH   =

#----------------------------------------------------------------------------
#	Compiler and linker commands
#----------------------------------------------------------------------------

FCOMP   = mpif90
CCOMP   = mpicc
CPPCOMP = mpiCC	
LINK    = mpif90

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

# for speed
FFLAGS_OPT = -fast -r8 -c

# for debugging
FFLAGS_DEBUG = -g -c 

F90FLAGS =


# C flags for Libraries
CFLAGS_HDF5 = -I $(HDF5_PATH)/include

# C flags for optimization
CFLAGS_OPT = -c -fast -DTFLOPS

# C flags for debugging
CFLAGS_DEBUG = -c -DTFLOPS


#----------------------------------------------------------------------------
# Linker flags
#----------------------------------------------------------------------------

# Linker flags for optimization
LFLAGS_OPT = -o

# Linker flags for debugging
LFLAGS_DEBUG = -g -o

#----------------------------------------------------------------------------
# Library specific linking
#----------------------------------------------------------------------------
LIB_OPT   =
LIB_DEBUG =

LIB_HDF5  = -L$(HDF5_PATH)/lib -lhdf5 -L$(ZLIB_PATH)/lib -lz -lnoop_stubs 
LIB_HDF4  = -L$(HDF4_PATH)/lib -lmfhdf -ldf -lz -lmpi

LIB_NCMPI =
LIB_MPE   =

#----------------------------------------------------------------------------
#	Additional machine-dependent object files
#----------------------------------------------------------------------------

MACHOBJ = 

#----------------------------------------------------------------------------
#	Additional commands
#----------------------------------------------------------------------------
MV   = mv -f
AR   = ar -r
RM   = rm -f
CD   = cd
RL   = ranlib
ECHO = echo
