#	FLASH makefile definitions for ix86 Linux (Portland Group compiler)

#	2/25/99

#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

MPI_HOME = /usr/mpich-1.2.4-ifc

#----------------------------------------------------------------------------
# Compiler and linker commands
#----------------------------------------------------------------------------

FCOMP   = mpif90
CCOMP   = cc
CPPCOMP = CC
LINK    = mpif90

# pre-processor flags
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

FFLAGS_OPT   =  -c -fast -r8 -i4 -I$(MPI_HOME)/include
FFLAGS_DEBUG =  -c -g -r8   -i4 -I$(MPI_HOME)/include
FFLAGS_TEST  =  -c -r8  -i4  -I$(MPI_HOME)/include

CFLAGS       = -c -I$(MPI_HOME)/include

#----------------------------------------------------------------------------
# Linker flags
#----------------------------------------------------------------------------

#  Linker flags for optimization
LFLAGS_OPT =  -o

#  Linker flags for debugging
LFLAGS_DEBUG = -g -o

#  Linker flags for testing
LFLAGS_TEST = -o

#----------------------------------------------------------------------------
# Library specific linking
#----------------------------------------------------------------------------

LIB_TEST = -L$(MPI_HOME)/lib -L/home/siegela/4.1r3_linux/lib/ \
           -lfmpich -lmpich -ljpeg -lz

LIB_OPT  = -L$(MPI_HOME)/lib -L/home/siegela/4.1r3_linux/lib/ \
           -lfmpich -lmpich -ljpeg -lz

LIB_DEBUG = -L$(MPI_HOME)/lib -L/home/siegela/4.1r3_linux/lib/ \
           -lfmpich -lmpich -ljpeg -lz

LIB_HDF5 = -lhdf5

#----------------------------------------------------------------------------
# Additional machine-dependent object files
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










