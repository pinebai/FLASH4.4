#FLASH makefile definitions for ix86 Linux (Portland Group compiler)
#
# 
# Portland Group Fortran 7.2-5 32-bit target on x86 Linux -tp k8-32 
# hdf5-1.8.7 (gcc 4.4.4)
# mpich-1.2.1p1 (gcc 4.4.4)

#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

HDF5_PATH = /opt/hdf5/pgi/1.8.8
MPI_PATH  = /opt/mpich2/pgi/1.2.1p1
PG_PATH   = /opt/pgi/linux86/7.2

HYPRE_PATH = /opt/hypre/pgi/2.8.0b
ifeq ("$(USEOPENMP)", "1")
HYPRE_PATH=/opt/hypre/pgi/2.8.0b_omp
endif



PAPI_PATH  = 
PAPI_FLAGS = 

FISHPAK_PATH = 
 

ZLIB_PATH  =

NCMPI_PATH = /opt/netcdf/pgi/current
MPE_PATH   =

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#  We use the f90 compiler as the linker, so some C libraries may explicitly
#  need to be added into the link line.
#----------------------------------------------------------------------------

#FCOMP   =  $(PG_PATH)/bin/pgf90
#CCOMP   =  $(PG_PATH)/bin/pgcc
#CPPCOMP =  $(PG_PATH)/bin/pgCC
#LINK    =  $(PG_PATH)/bin/pgf90

FCOMP      = $(MPI_PATH)/bin/mpif90
CCOMP      = $(MPI_PATH)/bin/mpicc
CPPCOMP    = $(MPI_PATH)/bin/mpicxx
LINK       = $(MPI_PATH)/bin/mpif90


# Any additional machine specific libraries to be linked in

# pre-processor flag

PP         = -D

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


#A note about Portland Group:
#the flags -fastsse and -pc 64 are used to avoid an 80-bit floating point 
#promotion that is a holdover from the x87 math co-processor instruction set

OPENMP = -mp -Minfo=mp

FFLAGS_OPT   = -c -r8 -i4 -fastsse -Mnovect -pc 64 
FFLAGS_DEBUG = -c -r8 -i4 -g -pc 64 -Ktrap=divz -Mbounds
#Ktrap=divz -Mchkfpstk -Mchkptr -Mchkstk -Mbounds
FFLAGS_TEST  = ${FFLAGS_OPT} -Kieee
#FFLAGS_TEST  = -c -r8 -i4 -fast

FFLAGS_HYPRE = -I${HYPRE_PATH}/include
CFLAGS_HYPRE = -I${HYPRE_PATH}/include


CFLAGS_OPT   = -c -O2
CFLAGS_DEBUG = -c -g
CFLAGS_TEST  = ${CFLAGS_OPT} -Kieee
#CFLAGS_TEST  = -c 

FFLAGS_PAPI  = -I$(PAPI_PATH)/include
FFLAGS_MPI   = -I$(MPI_PATH)/include

# if we are using HDF5, we need to specify the path to the include files
CFLAGS_HDF5  = -I$(HDF5_PATH)/include -DH5_USE_16_API
CFLAGS_NCMPI = -I$(NCMPI_PATH)/include
CFLAGS_MPI   = -I$(MPI_PATH)/include

#CFLAGS_SAMRAI = -I$(SAMRAI_PATH)/include

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT, 
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT   = -o
LFLAGS_DEBUG = -Wl,-noinhibit-exec -g -o
LFLAGS_TEST  = -Mprof=lines -o

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


LIB_HDF5    = -L$(HDF5_PATH)/lib -lhdf5 -lz  -lpthread
              
LIB_MPI     = -L$(MPI_PATH)/lib -lmpich
LIB_PNG     = -lpng -lz

LIB_OPT     = 
LIB_DEBUG   =
LIB_TEST    =

LIB_NCMPI   = -L$(NCMPI_PATH)/lib -lpnetcdf
LIB_MPE     =

LIB_HYPRE   = -L${HYPRE_PATH}/lib -lHYPRE

#LIB_SAMRAI  = -L$(SAMRAI_PATH)/lib -lSAMRAI
#LIB_PAPI    = $(PAPI_PATH)/lib/libpapi.a $(PAPI_PATH)/lib/_fixunssfdi.o

#LIB_FISHPAK = -L$(FISHPAK_PATH)/lib -lfishpak
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

MV    = mv -f
AR    = ar -r
RM    = rm -f
CD    = cd
RL    = ranlib
ECHO  = echo
AWK   = awk
CAT   = cat

#----------------------------------------------------------------------------
# Specific compilation.
#  Workarounds for some compiler bugs.
#----------------------------------------------------------------------------


#Use the totalview memory debugger (we overwrite certain Makefile.h variables):

#Specify MEMORY_DEBUG=YES in order to link the totalview memory debugging library.
#When you load totalview, click Help->About Totalview, and check the version matches 
#the $(TVIEW) variable.
#MEMORY_DEBUG=NO
ifeq ($(MEMORY_DEBUG), YES)

#TVIEW = /usr/local/toolworks/totalview.8.4.1-5/linux-x86

#Totalview crashes when we link against the static pthread library.
#Work around is to use shared object pthread library (reported on totalview website).

LIB_HDF5    = -L$(HDF5_PATH)/lib -lhdf5 -lz -Bdynamic -lpthread

endif


#ifeq ($(FLASHBINARY),true)
#FFLAGS_WO_FASTSSE = $(patsubst -fastsse,-fast,$(FFLAGS))
# Compile the following files with the same flags as others (which depend on
# whether -opt, -debug, or -test was in effect for setup), except that
# -fastsse is replaced by -fast. This allows compilation ofthese files
# on the specific compiler version used for testing on zingiber.uchicago.edu:
#        pgf90 6.0-4 32-bit target on x86 Linux
# Without this workaround, the compiler generates invalid assembler code
# when the preprocessor symbol NFLUXES has the value 0.
#Grid_putFluxData.o Grid_getFluxData.o: %.o : %.F90
#	$(FCOMP) $(FFLAGS_WO_FASTSSE) $(F90FLAGS) $(FDEFINES) $<
#endif
