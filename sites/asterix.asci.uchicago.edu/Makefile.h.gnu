# FLASH makefile definitions for x86-64 Linux (GNU compilers)
#----------------------------------------------------------------------------
# Set the HDF5/MPI library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

MPI_PATH   = /opt/openmpi-1.8.6_gcc
HDF4_PATH  =
HDF5_PATH  = /opt/hdf5/gcc/1.8.16
HYPRE_PATH = /opt/hypre/gcc/2.9.0b

ZLIB_PATH  =

PAPI_PATH  =
PAPI_FLAGS =

FISHPAK_PATH =

NCMPI_PATH = /usr/local/netcdf
MPE_PATH   =

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#
#----------------------------------------------------------------------------

FCOMP   = ${MPI_PATH}/bin/mpif90
CCOMP   = ${MPI_PATH}/bin/mpicc
CPPCOMP = ${MPI_PATH}/bin/mpicxx
LINK    = ${MPI_PATH}/bin/mpif90

# pre-processor flag
PP      = -D

#----------------------------------------------------------------------------
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

OPENMP = -fopenmp

FFLAGS_OPT = -ggdb -c -O2 -fdefault-real-8 -fdefault-double-8 \
-Wuninitialized

#I explictly add -O0 because I found that compiling source files without
#an optimization flag generates the same object code as compiling source
#files with -O2.  The -O0 is required so that gdb no longer shows
#"<value optimized out>" for certain function arguments.

FFLAGS_DEBUG = -ggdb -c -O0 -fdefault-real-8 -fdefault-double-8 \
-pedantic -Wall -Waliasing \
-Wsurprising -Wconversion -Wunderflow \
-ffpe-trap=invalid,zero,overflow -fbacktrace -fbounds-check \
-fimplicit-none -fstack-protector-all

FFLAGS_TEST = -ggdb -c -O0 -fdefault-real-8 -fdefault-double-8

FFLAGS_HYPRE = -I${HYPRE_PATH}/include
CFLAGS_HYPRE = -I${HYPRE_PATH}/include


F90FLAGS =


#The macro _FORTIFY_SOURCE adds some lightweight checks for buffer
#overflows at both compile time and run time (only active at -O1 or higher)
#http://gcc.gnu.org/ml/gcc-patches/2004-09/msg02055.html
CFLAGS_OPT = -ggdb -c -O2 -Wuninitialized -D_FORTIFY_SOURCE=2

CFLAGS_DEBUG = -ggdb -c -O0 -Wno-div-by-zero -Wundef \
-Wconversion -Wstrict-prototypes -Wunreachable-code \
-pedantic -Wall -Winit-self -ftree-vrp -Wfloat-equal \
-Wunsafe-loop-optimizations -Wpadded -fstack-protector-all

CFLAGS_TEST = -ggdb -O0 -c


# if we are using HDF5, we need to specify the path to the include files
CFLAGS_HDF5 = -I${HDF5_PATH}/include -DH5_USE_16_API
CFLAGS_NCMPI = -I${NCMPI_PATH}/include
$?$ needed? $?$ CFLAGS_MPI   = -I$(MPI_PATH)/include

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT,
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT   = -ggdb -O2 -o
LFLAGS_DEBUG = -ggdb -O0 -o
LFLAGS_TEST  = -ggdb -O0 -o

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

LIB_OPT   =
LIB_DEBUG =
LIB_TEST  =

LIB_HDF4  =
LIB_HDF5 = -L $(HDF5_PATH)/lib -lhdf5_fortran -lhdf5 -lz

LIB_PAPI  =
LIB_MATH  =

LIB_MPI   =
LIB_NCMPI = -L ${NCMPI_PATH}/lib -lpnetcdf
LIB_MPE   =

#LIB_HYPRE = -L${HYPRE_PATH}/lib -lHYPRE -llapack -lblas
LIB_HYPRE = -L${HYPRE_PATH}/lib -lHYPRE -llapack -lblas -static
LIB_LIBNBC = -L/home/cdaley/software/libNBC/1.1.1/mpich-1.4.1p1_gnu/lib -lnbc

LIB_STDCXX = -lstdc++

# Uncomment (and change) the following line to use an external library file for double precision SPECFUN code.
#LIB_SPECFUN_DP = -lspecfun_dp

# Uncomment the following line to use electic fence memory debugger.
# Need the following environmental variable (see env.sh):
# export EF_ALLOW_MALLOC_0=1
#CONFIG_LIB = -L/usr/lib64 -lefence


LIB_FISHPAK =

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

FFLAGS_WO_WARNALL = $(patsubst -warn all,,$(FFLAGS))

#Turn off compiler error messages for paramesh files that use wrapper
#functions such as MPI_int_SSEND.
amr_migrate_tree_data.o : %.o : %.F90
	$(FCOMP) $(FFLAGS_WO_WARNALL) $(F90FLAGS) $(FDEFINES)	$<
mpi_amr_test_neigh_values.o : %.o : %.F90
	$(FCOMP) $(FFLAGS_WO_WARNALL) $(F90FLAGS) $(FDEFINES)	$<
mpi_amr_checkpoint_default.o : %.o : %.F90
	$(FCOMP) $(FFLAGS_WO_WARNALL) $(F90FLAGS) $(FDEFINES)	$<
mpi_amr_morton.o : %.o : %.F90
	$(FCOMP) $(FFLAGS_WO_WARNALL) $(F90FLAGS) $(FDEFINES)	$<

#Fortran 77 source lines exceed 72 characters
umap.o : %.o : %.F
	$(FCOMP) $(FFLAGS_WO_WARNALL) $(FDEFINES)	$<
fftsg.o : %.o : %.f
	$(FCOMP) $(FFLAGS_WO_WARNALL) $(FDEFINES)	$<
fftsg3d.o : %.o : %.f
	$(FCOMP) $(FFLAGS_WO_WARNALL) $(FDEFINES)	$<

#Files mix and match assumed shape arrays, assumed size arrays
#and scalars in function calls.  This is fine but it is viewed as
#a problem when using strict type checking compiler options.
fftpack.o : %.o : %.f90
	$(FCOMP) $(FFLAGS_WO_WARNALL) $(FDEFINES)	$<
gr_pfftDcftForward.o : %.o : %.F90
	$(FCOMP) $(FFLAGS_WO_WARNALL) $(FDEFINES)	$<
gr_pfftDcftInverse.o : %.o : %.F90
	$(FCOMP) $(FFLAGS_WO_WARNALL) $(FDEFINES)	$<

#ifort version 12.1.0 hangs during compilation of hy_ppm_sweep.F90
#unless we use -O0
hy_ppm_sweep.o : %.o : %.F90
	$(FCOMP) $(FFLAGS) -O0 $(F90FLAGS) $(FDEFINES)	$<

#ifort version 12.1.0 generates bad code for Grid_advanceDiffusion.F90
#from Grid/GridSolvers/HYPRE/Grid_advanceDiffusion.F90 when we use the
#-openmp option: this is very strange because this file contains no openmp.
FFLAGS_WO_OPENMP = $(patsubst $(OPENMP),,$(FFLAGS))
Grid_advanceDiffusion.o : %.o : %.F90
	$(FCOMP) $(FFLAGS_WO_OPENMP) $(F90FLAGS) $(FDEFINES)	$<

# The following GNU make special prevents that apparent dependencies
# on the file iso_c_binding.mod, which does usually not actually exist
# in the object directory but refers to the ISO_C_BINDING module known
# to FORTRAN compilers internally, trigger unnecessary recompilation
# of files that refer to the ISO_C_BINDING module.
.SECONDARY: iso_c_binding.mod

endif
