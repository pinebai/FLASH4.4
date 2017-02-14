# FLASH makefile definitions for SGI/IRIX platform

#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------
HDF4_PATH  = /opt/pkgs/HDF/4.1r2_irix64v6.4-n32
HDF5_PATH  = /opt/pkgs/hdf5-1.4.5-64

ZLIB_PATH  = /opt/pkgs/zlib-1.1.3

MPE_PATH   = /home/tonychan/mpe_work/install_SGIMPI_64
NCMPI_PATH = /opt/pkgs/parallel-netcdf64

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#  We use the f90 compiler as the linker, so some C libraries may explicitly
#  need to be added into the link line.
#----------------------------------------------------------------------------

FCOMP   = f90 -64
CCOMP   = cc -64
CPPCOMP = CC -64
LINK    = f90 -64

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

FFLAGS_OPT   = -c -Ofast=ip27 -OPT:Olimit=0:IEEE_arithmetic=3:roundoff=3 -IPA \
               -r8 -d8 -i4 -cpp -r10000 -LNO
FFLAGS_DEBUG = -c -DEBUG:subscript_check=ON:verbose_runtime=ON -r8 -d8 -i4 \
               -cpp -g 
FFLAGS_TEST  = -c -r8 -d8 -i4 -cpp -O2

F90FLAGS = 

CFLAGS_OPT = -IPA -Ofast=ip27 -c
CFLAGS_DEBUG = -g -c
CFLAGS_TEST = -c -O2

CFLAGS_HDF5  = -I $(HDF5_PATH)/include

CFLAGS_NCMPI = -I $(NCMPI_PATH)/include
CFLAGS_VISTOOLS = -I /scratch2/caceres/python/include/python2.1

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT, 
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT   = -r8 -d8 -i4 -IPA -o
LFLAGS_DEBUG = -r8 -d8 -i4 -g -o
LFLAGS_TEST  = -r8 -d8 -i4 -o


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
LIB_OPT      = -L/usr/lib64 -lmpi -lfastm
LIB_DEBUG    = -L/usr/lib64 -lmpi -lfastm
LIB_TEST     = -L/usr/lib64 -lmpi 

LIB_HDF4     = -L$(HDF4_PATH)/lib -lmfhdf -ldf -lz

LIB_HDF5     = -B static -L $(HDF5_PATH)/lib -lhdf5 \
               -B dynamic -L $(ZLIB_PATH)/lib -lz

LIB_VISTOOLS = -L /scratch2/caceres/python/src -lpython2.1 \
               -lpthread
LIB_MATH     = -lcomplib.sgimath

LIB_NCMPI    = -B static -L $(NCMPI_PATH)/lib -lnetcdf

LIB_MPE      = -L$(MPE_PATH)/lib -lmpe_f2cmpi -llmpe -lmpe


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
