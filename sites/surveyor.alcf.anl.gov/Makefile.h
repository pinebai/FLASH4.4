#-------------------------------------------------------------------
# FLASH makefile definitions for BG/P (Surveyor)
#-------------------------------------------------------------------

#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

HDF5_PATH = /soft/apps/hdf5-1.8.0

ZLIB_PATH =

HPM_PATH =
PNG_PATH = 

NCMPI_PATH = /soft/apps/parallel-netcdf-1.0.3-xl
MPI_PATH = 

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#  We use the f90 compiler as the linker, so some C libraries may explicitly
#  need to be added into the link line.
#----------------------------------------------------------------------------

ifdef PDTDIR
OPT=-optPDBFile=merged.pdb -optTauSelectFile="select.tau" -optReset="" -optVerbose
else
#LIB_MPI = -L/$(MPI_PATH)/lib -lmpich.rts -lmsglayer.rts -ldevices.rts -lrts.rts -ldevices.rts
LIB_MPI = 
endif
#FCOMP   = mpif90.ibm
#CCOMP   = mpicc.ibm -I${MPI_PATH}/include
#CPPCOMP = mpif90.ibm -I${MPI_PATH}/include
#LINK    = mpif90.ibm

#Rob has built his own MPI library that includes an MPI-IO fix needed for HDF5
#collective IO.  Also this version does not leak memory when communicators 
#are destroyed meaning collective network does not need to be disabled.
FCOMP    = ~robl/soft/dcmf-2009078/bin/mpixlf90_r
CCOMP    = ~robl/soft/dcmf-2009078/bin/mpixlc_r
CPPCOMP  = ~robl/soft/dcmf-2009078/bin/mpixlcxx_r
LINK     = ~robl/soft/dcmf-2009078/bin/mpixlf90_r


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

FFLAGS_OPT   = -g -O2 -qintsize=4 -qrealsize=8 -qfixed -qnosave -c \
               -qinline -qmaxmem=16384 \
               -qsuffix=cpp=F -qarch=450 -qtune=auto -qcache=auto 
               
FFLAGS_TEST  = -g -O -qintsize=4 -qrealsize=8 -qfixed -qnosave -c \
               -qsuffix=cpp=F -qarch=450 -qtune=auto -qcache=auto -qmaxmem=16384
              
FFLAGS_DEBUG = -O -g -qintsize=4 -qrealsize=8 -qfixed -qnosave -c \
               -qarch=450 -qmaxmem=16384 -qcheck

F90FLAGS     = -qsuffix=f=F90:cpp=F90 -qfree=f90
f90FLAGS     = -qsuffix=f=f90:cpp=F90 -qfree=f90

# if we are using HDF5, we need to specify the path to the include files

CFLAGS_OPT   = -g -O2 -DIBM -DNOUNDERSCORE -c \
               -qarch=450 -qtune=auto -qcache=auto -qmaxmem=16384 -D_FILE_OFFSET_BITS=64
CFLAGS_TEST  = -g -O -DIBM -DNOUNDERSCORE -c \
               -qarch=450 -qtune=auto -qcache=auto -qmaxmem=16384 -D_FILE_OFFSET_BITS=64
CFLAGS_DEBUG = -g  -DIBM -DNOUNDERSCORE -c -qarch=450 -qmaxmem=16384 -D_FILE_OFFSET_BITS=64 -qcheck=all -qdbxextra

CFLAGS_HDF5  = -I${HDF5_PATH}/include -DNOUNDERSCORE -DH5_USE_16_API
CFLAGS_NCMPI = -I$(NCMPI_PATH)/include

MDEFS = -WF,

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


#----------------------------------------------------------------------------
# Specific compilation.
#  Can't compile certain files with array bounds checking because: 
#  "Zero-sized arrays must not be subscripted"
#---------------------------------------------------------------------------- 
Grid_getBlkData.o : Grid_getBlkData.F90
	${FCOMP} ${FFLAGS_TEST} ${F90FLAGS} ${FDEFINES} $<

Grid_getPlaneData.o : Grid_getPlaneData.F90
	${FCOMP} ${FFLAGS_TEST} ${F90FLAGS} ${FDEFINES} $<

Grid_getPointData.o : Grid_getPointData.F90
	${FCOMP} ${FFLAGS_TEST} ${F90FLAGS} ${FDEFINES} $<

Grid_getRowData.o : Grid_getRowData.F90
	${FCOMP} ${FFLAGS_TEST} ${F90FLAGS} ${FDEFINES} $<

Grid_putBlkData.o : Grid_putBlkData.F90
	${FCOMP} ${FFLAGS_TEST} ${F90FLAGS} ${FDEFINES} $<

Grid_putPlaneData.o : Grid_putPlaneData.F90
	${FCOMP} ${FFLAGS_TEST} ${F90FLAGS} ${FDEFINES} $<

Grid_putPointData.o : Grid_putPointData.F90
	${FCOMP} ${FFLAGS_TEST} ${F90FLAGS} ${FDEFINES} $<

Grid_putRowData.o : Grid_putRowData.F90
	${FCOMP} ${FFLAGS_TEST} ${F90FLAGS} ${FDEFINES} $<

IO_init.o : IO_init.F90
	${FCOMP} ${FFLAGS_TEST} ${F90FLAGS} ${FDEFINES} $<

io_writeData.o : io_writeData.F90
	${FCOMP} ${FFLAGS_TEST} ${F90FLAGS} ${FDEFINES} $<

cma_flatten.o : cma_flatten.F90
	${FCOMP} ${FFLAGS_TEST} ${F90FLAGS} ${FDEFINES} $<

hydro_1d.o : hydro_1d.F90
	${FCOMP} ${FFLAGS_TEST} ${F90FLAGS} ${FDEFINES} $<

hy_ppm_block.o : hy_ppm_block.F90
	${FCOMP} ${FFLAGS_TEST} ${F90FLAGS} ${FDEFINES} $<

hy_ppm_updateSoln.o : hy_ppm_updateSoln.F90
	${FCOMP} ${FFLAGS_TEST} ${F90FLAGS} ${FDEFINES} $<

intrfc.o : intrfc.F90
	${FCOMP} ${FFLAGS_TEST} ${F90FLAGS} ${FDEFINES} $<

rieman.o : rieman.F90
	${FCOMP} ${FFLAGS_TEST} ${F90FLAGS} ${FDEFINES} $<

states.o : states.F90
	${FCOMP} ${FFLAGS_TEST} ${F90FLAGS} ${FDEFINES} $<


# add_block_to_tree subroutine in local_tree.F90 is called with the same argument twice.  
# The arguments alias each other so we have to inform the compiler not to use unsafe 
# optimisations (-qalias).  Required at -O4 optimisation level for correctness.
local_tree_module.mod local_tree.mod local_tree.o : local_tree.F90
	${FCOMP} ${FFLAGS} -qalias=nostd ${F90FLAGS} ${FDEFINES} $<
