#       FLASH makefile definitions for TRU64/Alpha (Compaq compiler, MPICH)
#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

#HDF4_PATH = /usr/local/packages/HDF4.1r4
#HDF5_PATH = /usr/local/packages/HDF5-1.4.2

DPMTA_PATH = /users/zuhone/dpmta-3.1.3

ZLIB_PATH  = 

PAPI_PATH  = 
PAPI_FLAGS = 

MPI_PATH   = /usr/local/opt/Compaq_MPI_64_2.6_r6

MPE_PATH   = 
NCMPI_PATH =

LIB_PNG    = 

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#  We use the f90 compiler as the linker, so some C libraries may explicitly
#  need to be added into the link line.
#----------------------------------------------------------------------------

FCOMP   = f95
CCOMP   = cc
CPPCOMP = CC
LINK    = f95

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

FFLAGS_OPT      = -c -fast -cord -p -pg -inline speed -math_library fast -O5 -r8 -i4 -cpp \
                  -fpe1 -double_size 64 $(MPI_COMPILE_FLAGS) $(HDF_COMPILE_FLAGS)
#FFLAGS_DEBUG   = -c -g -pg -p1 -check bounds -fuse_xref -ladebug -warn argument_checking \
                  -warn declarations -warn unused  -r8 -i4 -cpp -fpe1 -double_size 64 \
                  $(MPI_COMPILE_FLAGS) $(HDF_COMPILE_FLAGS)
FFLAGS_DEBUG    = -c -g -pg -p1 -check bounds -fuse_xref -ladebug -warn argument_checking \
                  -warn declarations -warn unused  -r8 -i4 -cpp -fpe1 -double_size 64 \
                  $(MPI_COMPILE_FLAGS) $(HDF_COMPILE_FLAGS)
FFLAGS_TEST     = -c -fast -check bounds -fuse_xref -ladebug -warn argument_checking \
                  -warn declarations -warn unused  -r8 -i4 -cpp -fpe1 -double_size 64 \
                  $(MPI_COMPILE_FLAGS) $(HDF_COMPILE_FLAGS)

F90FLAGS        =


CFLAGS_OPT      = -I. -I/usr/include -fast -inline speed -math_library fast -O5 -fpe1 \
                  -c $(MPI_COMPILE_FLAGS) $(HDF_COMPILE_FLAGS)
CFLAGS_DEBUG    = -I. -I/usr/include -g -pg -pg1 -check bounds -fuse_xref -ladebug -fpe1 \
                  -DDEBUG_AMR -c $(MPI_COMPILE_FLAGS) $(HDF_COMPILE_FLAGS)
#CFLAGS_DEBUG   = -I. -I/usr/include -g -pg -pg1 -check bounds -fuse_xref -ladebug -fpe1 \
#                  -c $(MPI_COMPILE_FLAGS) $(HDF_COMPILE_FLAGS)
CFLAGS_TEST     = -I. -I/usr/include -fast -fuse_xref -fpe1 -c $(MPI_COMPILE_FLAGS) \
                  $(HDF_COMPILE_FLAGS)
 
#CFLAGS_HDF5    = -I$(HDF5_PATH)/include
CFLAGS_HDF5     = $(HDF_COMPILE_FLAGS)
CFLAGS_DPMTA    = -I$(DPMTA_PATH)/include

CFLAGS_NCMPI    =
CFLAGS_VISTOOLS = 

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT,
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT   = -r8 -double_size 64 -i4 -IPA $(MPI_LD_FLAGS) $(HDF_LD_FLAGS) -o
LFLAGS_DEBUG = -r8 -double_size 64 -i4 -g $(MPI_LD_FLAGS) $(HDF_LD_FLAGS) -o
LFLAGS_TEST  = -r8 -double_size 64 -i4 $(MPI_LD_FLAGS) $(HDF_LD_FLAGS) -o


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
LIB_HDF4  = -lmfhdf -ldf -lz
LIB_HDF5  = -lmpio -lhdf5 -lz
LIB_DPMTA = -L$(DPMTA_PATH)/lib -ldpmta

MPI_LIB   = 
LIB_MPI   =

LIB_MPE   = 
LIB_PAPI  = 
LIB_NCMPI = 

LIB_OPT   =  -lfmpi -lmpi -lelan -math_library fast
LIB_DEBUG =  -lfmpi -lmpi -lelan -math_library fast
LIB_TEST  =  -lfmpi -lmpi -lelan -math_library fast


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
