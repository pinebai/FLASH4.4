# This Makefile will use the macro definitions from the
# FLASH Makefile.h in the sites directory.  The site
# directory name is obtained from the hostname command.
#
# To build sfocu:
# > make
#
# Alternative ways to build sfocu if a simple make fails:
#
# 1. If the sites directory has a different name to the return
# value of hostname command then you must define SITE macro to be
# the sites directory of interest:
# > make SITE=code.uchicago.edu
#
#
# 2. If you do not have a Pnetcdf installation then you can still
# build sfocu using:
# > make NO_NCDF=True
#
#
# 3. If you only have a serial installation of HDF5 then you
# can still build sfocu using:
# > make NO_MPI=True
#
# You also need to 1). specify NO_NCDF=True because Pnetcdf depends
# on MPI and 2). make CCOMP the actual compiler and not an MPI
# compiler wrapper script:
# > make NO_MPI=True NO_NCDF=True CCOMP=gcc
#
#
# 4. If you want to build sfocu for execution on Blue Gene login nodes
# then you must make use of a login node installation of HDF5.
#
# make SITE=cetuslac1.i2b.alcf.anl.gov NO_MPI=True NO_NCDF=True CCOMP=xlc_r \
#      CFLAGS_HDF5='-I/veas-fs0/cdaley/software/hdf5-fen/1.8.8/xl/include -DH5_USE_16_API' \
#      LIB_HDF5='-L/veas-fs0/cdaley/software/hdf5-fen/1.8.8/xl/lib -lhdf5'

# Include the appropriate Makefile.h.
ifdef SITE
	SITEDIR=$(SITE)
else
	SITEDIR:=$(shell hostname)
endif
include ../../sites/$(SITEDIR)/Makefile.h


# Use macros from Makefile.h.
ifdef NO_NCDF
	CDEFINES_NCMPI=-DNO_NCDF
	_LIB_NCMPI=
else
	CDEFINES_NCMPI=$(CFLAGS_NCMPI)
	_LIB_NCMPI=$(LIB_NCMPI)
endif

ifdef NO_MPI
	CDEFINES_MPI=
	_LIB_MPI=
else
	CDEFINES_MPI=$(CFLAGS_MPI) -DNEED_MPI
	_LIB_MPI=$(LIB_MPI)
endif


CDEFINES_HDF5=$(CFLAGS_HDF5) -DNO_HDF4
CDEFINES=$(CDEFINES_HDF5) $(CDEFINES_NCMPI) $(CDEFINES_MPI)

CFLAGS=$(CFLAGS_OPT)
LIBS=$(LIB_HDF5) $(_LIB_NCMPI) $(_LIB_MPI) -lm
LFLAGS=$(LFLAGS_OPT)


EXENAME = sfocu
OBJS = main.o sfocu.o test_reader.o flash_reader.o namecmp.o \
       flash_reader_hdf5.o flash_reader_ncdf.o flash_reader_chombo.o sameblock.o options.o

.SUFFIXES: .c .o
.c.o :
	$(CCOMP) $(CFLAGS) $(CDEFINES) $*.c

default: $(EXENAME)

$(EXENAME): $(OBJS)
	$(CCOMP) $(OBJS) $(LIBS) $(LFLAGS) $(EXENAME)

clean:
	rm -f $(OBJS) $(EXENAME)
