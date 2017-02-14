# This Makefile will use the macro definitions from the
# FLASH Makefile.h in the sites directory.  The site
# directory name is obtained from the hostname command.
#
# To build cmpExpIO:
# > make
#
# Alternative ways to build cmpExpIO if a simple make fails:
#
# 1. If the sites directory has a different name to the return
# value of hostname command then you must define SITE macro to be
# the sites directory of interest:
# > make SITE=code.uchicago.edu
#
# 2. If you do not have a Pnetcdf installation then you can still
# build cmpExpIO using:
# > make NO_NCDF=True
#

# Include Makefile.h.
ifdef SITE
	SITEDIR=$(SITE)
else
	SITEDIR:=$(shell hostname)
endif
include ../../sites/$(SITEDIR)/Makefile.h

# Use macros from Makefile.h.
ifdef NO_NCDF
	_CFLAGS_NCMPI=-DNO_NCDF
	_LIB_NCMPI=
	CPPSRC_NCDF=
else
	_CFLAGS_NCMPI=$(CFLAGS_NCMPI)
	_LIB_NCMPI=$(LIB_NCMPI)
	CPPSRC_NCDF=PNETCDF_File.cpp
endif


ifdef NO_HDF5
	_CFLAGS_HDF5=-DNO_HDF5
	_LIB_HDF5=
	CPPSRC_HDF5=
else
	_CFLAGS_HDF5=$(CFLAGS_HDF5)
	_LIB_HDF5=$(LIB_HDF5)
	CPPSRC_HDF5=HDF5_File.cpp
endif

CPPSRC	= WorkCoordinator.cpp \
	FlashFileFactory.cpp \
	DataBufferComparer.cpp \
	main.cpp \
	$(CPPSRC_NCDF) \
	$(CPPSRC_HDF5)

ifdef DEBUG
CFLAGS_OPT += -DDEBUG
endif

CFLAGS=$(CFLAGS_OPT) $(_CFLAGS_HDF5) $(_CFLAGS_NCMPI) $(CFLAGS_MPI) \
	-DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX
LIBS=$(_LIB_HDF5) $(_LIB_NCMPI) $(LIB_MPI)
LFLAGS=$(LFLAGS_OPT)


CPPOBJ	= $(CPPSRC:.cpp=.o)
EXE     = cmpExpIO

all:    $(CPPOBJ)
	$(CPPCOMP) ${LFLAGS} $(EXE) $(CPPOBJ) $(LIBS)

%.o:    %.cpp
	$(CPPCOMP) $< $(CFLAGS) -o $@

clean:
	rm -f $(CPPOBJ) $(EXE) core*
