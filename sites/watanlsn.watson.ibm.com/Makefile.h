#-------------------------------------------------------------------
# FLASH makefile definitions for BG/Q (EAS)
# Running make will use the raw compilers.
# We can use the MPI wrapper scripts with make USE_MPI_WRAPPERS=True
#-------------------------------------------------------------------

MPI_PATH = /bgsys/drivers/ppcfloor/comm/xl
MASS = -lmass

#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

HDF5_PATH = /home/xamorozov/hdf5-1.8.7
NCMPI_PATH = /home/xacdaley/software/parallel-netcdf/1.2.0/xl
ZLIB_PATH =
PNG_PATH =

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#  We use the f90 compiler as the linker, so some C libraries may explicitly
#  need to be added into the link line.
#----------------------------------------------------------------------------

ifdef USE_MPI_WRAPPERS
FCOMP    = $(MPI_PATH)/bin/mpixlf90_r
CCOMP    = $(MPI_PATH)/bin/mpixlc_r
CPPCOMP  = $(MPI_PATH)/bin/mpixlcxx_r
LINK     = $(MPI_PATH)/bin/mpixlf90_r
else
MPI_INC = -I$(MPI_PATH)/include
FCOMP    = /opt/ibmcmp/xlf/bg/14.1/bin/bgxlf90_r $(MPI_INC)
CCOMP    = /opt/ibmcmp/vacpp/bg/12.1/bin/bgxlc_r $(MPI_INC)
CPPCOMP  = /opt/ibmcmp/vacpp/bg/12.1/bin/bgxlc++_r $(MPI_INC)
LINK     = /opt/ibmcmp/xlf/bg/14.1/bin/bgxlf90_r
endif

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

OPENMP = -qsmp=omp
FFLAGS_OPT   = -g -O3 -qnohot -qrealsize=8 -qnosave -qfixed -c -qthreaded
FFLAGS_TEST  = $(FFLAGS_OPT)
FFLAGS_DEBUG = $(FFLAGS_OPT)

F90FLAGS     = -qsuffix=f=F90:cpp=F90 -qfree=f90
f90FLAGS     = -qsuffix=f=f90:cpp=F90 -qfree=f90

# if we are using HDF5, we need to specify the path to the include files

CFLAGS_OPT   = -g -O3 -qnohot -DIBM -c -D_FILE_OFFSET_BITS=64
CFLAGS_TEST  = $(CFLAGS_OPT)
CFLAGS_DEBUG = $(CFLAGS_OPT) -qcheck=all -qdbxextra

CFLAGS_HDF5  = -I$(HDF5_PATH)/include -DH5_USE_16_API
CFLAGS_NCMPI = -I$(NCMPI_PATH)/include

MDEFS = -WF,

.SUFFIXES: .o .c .f .F .h .fh .F90 .f90

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT,
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------
LFLAGS_OPT   = -g -O3 -qnohot -qnosave -o
LFLAGS_TEST  = -pg $(LFLAGS_OPT)
LFLAGS_DEBUG = $(LFLAGS_OPT)

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

ifdef USE_MPI_WRAPPERS
LIB_MPI =
else
LIB_MPI = -L$(MPI_PATH)/lib -lmpich -lmpl -lopa -lmpl \
          -L/bgsys/drivers/ppcfloor/comm/sys/lib -lpami \
          -L/bgsys/drivers/ppcfloor/spi/lib -lSPI_cnk -lrt -lpthread
endif

#If we included the mpihpm Profiler unit then LIB_MPIHPM or LIBMPIHPM_SMP
#is used in the actual Makefile.  If not then LIB_MPI is used.
HPM_COUNTERS = /bgsys/drivers/ppcfloor/bgpm/lib/libbgpm.a \
               /bgsys/drivers/ppcfloor/spi/lib/libSPI_upci_cnk.a
LIB_MPIHPM = -L/home/xamorozov/HPM/lib -lmpihpm $(HPM_COUNTERS) $(LIB_MPI)
LIB_MPIHPM_SMP = -L/home/xamorozov/HPM/lib -lmpihpm_smp $(HPM_COUNTERS) \
                 -L/opt/ibmcmp/xlsmp/bg/3.1/bglib64 -lxlsmp $(LIB_MPI)



LIB_HDF5  = -L$(HDF5_PATH)/lib -lhdf5
LIB_NCMPI = -L$(NCMPI_PATH)/lib -lpnetcdf
LIB_MATH  =

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
