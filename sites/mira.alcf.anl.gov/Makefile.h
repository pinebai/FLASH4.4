#-------------------------------------------------------------------
# FLASH makefile definitions for BG/Q (MIRA)
#-------------------------------------------------------------------

MASS = -lmass

#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

#CD: Use my installation of HDF5 because the public installation at
#/soft/libraries has not been rebuilt since Feb 15 2012 and therefore
#may have issues with recent driver upgrades
HDF5_PATH = /projects/Omega-NIF_Exp/tzeferac/software/V1R2M0_efix32/hdf5/1.8.10/xl
HYPRE_PATH = /projects/Omega-NIF_Exp/tzeferac/software/V1R2M0_efix32/hypre/2.8.0b/xl
#HDF5_PATH = /soft/libraries/hdf5/1.8.10/cnk-xl/current/ 
#/soft/libraries/3rdparty/hdf5-1.8.8

NCMPI_PATH =
ZLIB_PATH =
PNG_PATH =

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#  We use the f90 compiler as the linker, so some C libraries may explicitly
#  need to be added into the link line.
#----------------------------------------------------------------------------

FCOMP    = mpixlf90_r
CCOMP    = mpixlc_r
CPPCOMP  = mpixlcxx_r
LINK     = mpixlf90_r
#Note the plugin option for hpctoolkit which brings in a new feature to
#find the location in a program where threads are idle.
#LINK     = /soft/perftools/hpctoolkit/bin/hpclink --plugin omp-ibmxl mpixlf90_r

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

OPENMP = -qsmp=omp:noauto

#It is better to use a quiet NAN instead of a signaling NAN because
#the signaling NAN triggers too many false positives.  There are many
#places in FLASH where we copy around uninitialized data in array
#copies like a = b.  In such a copy there can be some uninitialized
#elements but they never propogate into the end numeric solution.

#Note: do not link against libmass if you want to catch a floating
#point invalid operation in a log (e.g. log(-1)).
SIGNALING_NAN='7ff7ffff'
QUIET_NAN='ff'
NAN=${QUIET_NAN}

FFLAGS_OPT   = -g -O3 -qnohot -qrealsize=8 -qnosave -qfixed -c -qthreaded
FFLAGS_TEST  = -g -O3 -qnohot -qstrict=all -qrealsize=8 -qnosave -qfixed -c -qthreaded
FFLAGS_DEBUG = -g -qnoopt -qrealsize=8 -qnosave -qfixed -c -qthreaded \
               -qinitauto=${NAN} -qinitalloc=${NAN} \
               -qflttrap=enable:invalid:overflow:zerodivide -qsigtrap=xl__trcedump
FFLAGS_HYPRE = -I${HYPRE_PATH}/include

F90FLAGS     = -qsuffix=f=F90:cpp=F90 -qfree=f90
f90FLAGS     = -qsuffix=f=f90:cpp=F90 -qfree=f90

# if we are using HDF5, we need to specify the path to the include files

CFLAGS_OPT   = -g -O3 -qnohot -DIBM -c -D_FILE_OFFSET_BITS=64
CFLAGS_TEST  = $(CFLAGS_OPT) -qstrict=all
CFLAGS_DEBUG = -g -qnoopt -DIBM -c -D_FILE_OFFSET_BITS=64 \
               -qcheck=all -qdbxextra \
               -qinitauto=${NAN} -qflttrap=enable:invalid:overflow:zerodivide

CFLAGS_HDF5  = -I$(HDF5_PATH)/include -DH5_USE_16_API
CFLAGS_NCMPI = -I$(NCMPI_PATH)/include
CFLAGS_HYPRE = -I${HYPRE_PATH}/include

MDEFS = -WF,

.SUFFIXES: .o .c .f .F .h .fh .F90 .f90

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT,
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------
LFLAGS_OPT   = -g -O3 -qnohot -qnosave -o
LFLAGS_TEST  = -g -O3 -qnohot -qnosave -o
LFLAGS_DEBUG = -g -qnoopt -qnosave -o

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

LIB_MPI =

#If we included the mpihpm Profiler unit then LIB_MPIHPM or LIBMPIHPM_SMP
#is used in the actual Makefile.  If not then LIB_MPI is used.
#I haven't tested mpihpm on Mira yet.
HPM_COUNTERS = /bgsys/drivers/ppcfloor/bgpm/lib/libbgpm.a
#LIB_MPIHPM = -L/soft/perftools/hpctw -lmpihpm $(HPM_COUNTERS) $(LIB_MPI)
#LIB_MPIHPM_SMP = -L/soft/perftools/hpctw -lmpihpm_smp $(HPM_COUNTERS) $(LIB_MPI)

LIB_MPIHPM = /home/morozov/HPM/lib/libmpihpm.a $(HPM_COUNTERS) $(LIB_MPI)
LIB_MPIHPM_SMP = /home/morozov/HPM/lib/libmpihpm_smp.a $(HPM_COUNTERS) $(LIB_MPI)
#Use cprof to create a profile file.
#~morozov/bin/cprof -e -n flash4 vmon.out.0 > profile.0


LIB_HDF5  = -L$(HDF5_PATH)/lib -lhdf5
LIB_NCMPI = -L$(NCMPI_PATH)/lib -lpnetcdf
LIB_MATH  =
LIB_HYPRE = -L${HYPRE_PATH}/lib -lHYPRE

LIB_OPT   = $(MASS)
LIB_DEBUG =
LIB_TEST  = $(MASS)

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

ifeq ($(FLASHBINARY),true)

# add_block_to_tree subroutine in local_tree.F90 is called with the same argument twice.
# The arguments alias each other so we have to inform the compiler not to use unsafe
# optimisations (-qalias).  Required at -O4 optimisation level for correctness.
local_tree_module.mod local_tree.mod local_tree.o : local_tree.F90
	${FCOMP} ${FFLAGS} -qalias=nostd ${F90FLAGS} ${FDEFINES} $<

endif
