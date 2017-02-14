#-------------------------------------------------------------------
# FLASH makefile definitions for the llnl bgl
# 
#
# These modules set environment variables pointing to the
# proper library locations. Only python is essential.
#-------------------------------------------------------------------

#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

HDF5_PATH = /usr/local/tools/hdf5/hdf5-1.6.4/parallel
HDF5_INCLUDE = -I/usr/local/tools/hdf5/hdf5-1.6.4/parallel/include

NCMPI_PATH = /usr/local/tools/parallel-netcdf/parallel-netcdf-1.0.0

ZLIB_PATH =

HPM_PATH =
PNG_PATH = /usr

MPI_PATH = #/bgl/BlueLight/ppcfloor/bglsys

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#  We use the f90 compiler as the linker, so some C libraries may explicitly
#  need to be added into the link line.
#----------------------------------------------------------------------------

ifdef PDTDIR
OPT=-optPDBFile=merged.pdb -optTauSelectFile="select.tau" -optReset="" -optVerbose
else
LIB_MPI = #-L$(MPI_PATH)/lib -lmpich.rts -lmsglayer.rts -ldevices.rts -lrts.rts -ldevices.rts
endif
FCOMP   = $(TAU_COMPILER) $(OPT) mpxlf90 -I/bgl/BlueLight/ppcfloor/bglsys/include
CCOMP   = $(TAU_COMPILER) mpxlc -I/bgl/BlueLight/ppcfloor/bglsys/include
CPPCOMP = $(TAU_COMPILER) mpxlC -I/bgl/BlueLight/ppcfloor/bglsys/include
LINK    = $(TAU_COMPILER) mpxlf90  

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

FFLAGS_OPT   = -O5 -qintsize=4 -qrealsize=8 -qfixed -qnosave -c \
               -qinline -qmaxmem=-1 -qxlf90=autodealloc \
               -qsuffix=cpp=F -qstrict -qlistopt -qtune=auto -qcache=auto 
               
FFLAGS_TEST  = -O -qintsize=4 -qrealsize=8 -qfixed -qnosave -c \
               -qsuffix=cpp=F -qarch=440 -qtune=auto -qcache=auto -qmaxmem=16384
              
FFLAGS_DEBUG = -g -qintsize=4 -qrealsize=8 -qfixed -qnosave -c \
               -qarch=440 -qmaxmem=16384

F90FLAGS     = -qsuffix=f=F90:cpp=F90 -qfree=f90
f90FLAGS     = -qsuffix=f=f90:cpp=F90 -qfree=f90

# if we are using HDF5, we need to specify the path to the include files
CFLAGS_HDF5  = -I${HDF5_PATH}/include
CFLAGS_NCMPI = -I$(NCMPI_PATH)/include

CFLAGS_OPT   = -O5 -DIBM -DNOUNDERSCORE -c -qlanglvl=extended \
               -qtune=auto -qcache=auto -qmaxmem=-1 -qstrict
CFLAGS_TEST  = -O -DIBM -DNOUNDERSCORE -c -qlanglvl=extended\
                -qtune=auto -qcache=auto -qmaxmem=16384
CFLAGS_DEBUG = -g  -DIBM -DNOUNDERSCORE -c -qlanglvl=extended\
               -qmaxmem=16384

MDEFS = -WF,

.SUFFIXES: .o .c .f .F .h .fh .F90 .f90

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT, 
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT   = -blpdata -qipa -o  # something to do with large page size -> better performance
LFLAGS_TEST  = -blpdata -o
LFLAGS_DEBUG = -blpdata -g -o

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

LIB_HDF5  = -L$(HDF5_PATH)/lib -lhdf5
LIB_NCMPI = -L$(NCMPI_PATH)/lib -lpnetcdf




LIB_OPT   = 
LIB_DEBUG = -L/usr/local/totalview/rs6000/lib/m/up/ -L/usr/local/totalview/rs6000/lib \
/usr/local/totalview/rs6000/lib/aix_malloctype64_5.o

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
