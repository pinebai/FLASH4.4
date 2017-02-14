#-------------------------------------------------------------------
# FLASH makefile definitions for the NERSC IBM SP2 (seabor) : 64-bit
# 
#
# NOTE: seaborg makes a lot of use of Modules, so it is 
# important to make sure that the following modules are
# loaded in order to compile FLASH:
#
#  python
#  hdf5_par_64/1.4.5-post2
#
# Optional modules:
#
#  GNU
#  MASS
#
# These modules set environment variables pointing to the
# proper library locations. Only python is essential.
#-------------------------------------------------------------------

#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

HDF5_PATH = /soft/tools/hdf5-1.6.2-rev2

ZLIB_PATH =

HPM_PATH =
PNG_PATH = 

MPI_PATH = /bgl/BlueLight/ppcfloor/bglsys/

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#  We use the f90 compiler as the linker, so some C libraries may explicitly
#  need to be added into the link line.
#----------------------------------------------------------------------------

ifdef PDTDIR
OPT=-optPDBFile=merged.pdb -optTauSelectFile="select.tau" -optReset="" -optVerbose
else
LIB_MPI = -L/$(MPI_PATH)/lib -lmpich.rts -lmsglayer.rts -ldevices.rts -lrts.rts -ldevices.rts
endif
FCOMP   = $(TAU_COMPILER) $(OPT) blrts_xlf90 -I/bgl/BlueLight/ppcfloor/bglsys/include
CCOMP   = $(TAU_COMPILER) blrts_xlc -I/bgl/BlueLight/ppcfloor/bglsys/include
CPPCOMP = $(TAU_COMPILER) blrts_xlC -I/bgl/BlueLight/ppcfloor/bglsys/include
LINK    = $(TAU_COMPILER) blrts_xlf90  
#CONFIG_LIB = -L/bgl/BlueLight/ppcfloor/bglsys/lib -lmpich.rts -lmsglayer.rts -ldevices.rts -lrts.rts -ldevices.rts
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

FFLAGS_OPT   = -O2 -qintsize=4 -qrealsize=8 -qfixed -qnosave -c \
               -qinline -qmaxmem=16384 \
               -qsuffix=cpp=F -qarch=440 -qtune=auto -qcache=auto 
               
FFLAGS_TEST  = -O2 -qintsize=4 -qrealsize=8 -qfixed -qnosave -c \
               -qsuffix=cpp=F -qarch=440 -qtune=auto -qcache=auto -qmaxmem=16384
              
FFLAGS_DEBUG = -O2 -g -qintsize=4 -qrealsize=8 -qfixed -qnosave -c \
               -qarch=440 -qmaxmem=16384

F90FLAGS     = -qsuffix=f=F90:cpp=F90 -qfree=f90
f90FLAGS     = -qsuffix=f=f90:cpp=F90 -qfree=f90

# if we are using HDF5, we need to specify the path to the include files
#CFLAGS_HDF5  = -I${HDF5_PATH}/include -I${PNG_PATH}/include
CFLAGS_HDF5  = -I${HDF5_PATH}/include -DNOUNDERSCORE

CFLAGS_OPT   = -O -DIBM -DNOUNDERSCORE -c \
               -qarch=440 -qtune=auto -qcache=auto -qmaxmem=16384 -D_FILE_OFFSET_BITS=64
CFLAGS_TEST  = -O2 -DIBM -DNOUNDERSCORE -c \
               -qarch=440 -qtune=auto -qcache=auto -qmaxmem=16384 -D_FILE_OFFSET_BITS=64
CFLAGS_DEBUG = -g  -DIBM -DNOUNDERSCORE -c \
               -qarch=440 -qmaxmem=16384 -D_FILE_OFFSET_BITS=64

MDEFS = -WF,

.SUFFIXES: .o .c .f .F .h .fh .F90 .f90

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT, 
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT   = -o
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


LIB_HDF5  = -L${HDF5_PATH}/lib -lhdf5

LIB_MATH  = 

LIB_OPT   = ${MASS}
LIB_DEBUG =
LIB_TEST  =

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
