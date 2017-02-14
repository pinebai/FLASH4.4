#-------------------------------------------------------------------
# FLASH makefile definitions for ALCF (argonne) production BG/P (intrepid)
#  started from seaborg (NERSC) makefile
#
#  This is for the GNU compiler.
# 
#
#-------------------------------------------------------------------

#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

HDF5_PATH = /soft/apps/hdf5-1.6.6

ZLIB_PATH =

HPM_PATH =
PNG_PATH = 

NCMPI_PATH = /bgusr/robl/soft/pnetcdf-20071019
MPI_PATH = 

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#  We use the f90 compiler as the linker, so some C libraries may explicitly
#  need to be added into the link line.
#----------------------------------------------------------------------------

LIB_MPI = 

FCOMP   = mpif90
CCOMP   = mpicc
CPPCOMP = mpicxx
LINK    = mpif90

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

FFLAGS_OPT = -c -O2 -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none
FFLAGS_DEBUG = -c -g -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -ffpe-trap=invalid,zero,overflow -fbounds-check
FFLAGS_TEST = -c -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none

CFLAGS_OPT = -O2  -c -DIBM -DNOUNDERSCORE
CFLAGS_DEBUG = -c -g -fbounds-check -DIBM -DNOUNDERSCORE
CFLAGS_TEST = -c -DIBM -DNOUNDERSCORE

CFLAGS_HDF5  = -I${HDF5_PATH}/include -DNOUNDERSCORE
CFLAGS_NCMPI = -I$(NCMPI_PATH)/include


.SUFFIXES: .o .c .f .F .h .fh .F90 .f90

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT, 
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT   = -O2 -o 
LFLAGS_TEST  = -o
LFLAGS_DEBUG = -g -o

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


LIB_MPI   = 
LIB_HDF5  = -L${HDF5_PATH}/lib -lhdf5
LIB_NCMPI = -L$(NCMPI_PATH)/lib -lpnetcdf

LIB_MATH  = 

LIB_OPT   = 
LIB_DEBUG =
LIB_TEST  =

LIB_STDCXX = -lstdc++

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

# full optimization does not work on runtime_parameters.F90
#runtime_parameters.o : runtime_parameters.F90
#	${FCOMP} ${FFLAGS_TEST}  $(F90FLAGS) $<

#amr_%.o : amr_%.F90
#	${FCOMP} ${FFLAGS_OPT_NOIPA} ${F90FLAGS} ${FDEFINES} $<
