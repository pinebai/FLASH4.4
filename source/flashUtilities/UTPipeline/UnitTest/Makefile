# To build the unit test, type:
#
# This Makefile uses macro definitions from your FLASH site Makefile.h,
# where, the site name is obtained from the hostname command.  The
# SITE macro is available if you wish to use different site:
# > make SITE=code.uchicago.edu ...
#
# The HEADER macro is available if you wish to use a FLASH Makefile.h
# with a non-default name:
# > make HEADER=Makefile.h.gnu ...
#
# The MPI3 macro specifies you wish to use MPI-3.
# > make MPI3=1 ...
# 
# The LIBNBC macro specifies you wish to use LibNBC for MPI-3.
# > make LIBNBC=1 ...
#
# Full example:
# > make LIBNBC=1 HEADER=Makefile.h.gnu


# Include the appropriate Makefile.h.
ifdef SITE
	SITEDIR=$(SITE)
else
	SITEDIR:=$(shell hostname -f)
endif

ifdef HEADER
	MAKEFILE_HEADER=$(HEADER)
else
	MAKEFILE_HEADER=Makefile.h
endif


include ../../../../sites/$(SITEDIR)/$(MAKEFILE_HEADER)


ifdef MPI3
	FFLAGS_MPI3 = $(MDEFS)-DFLASH_MPI3
endif

ifdef LIBNBC
	FFLAGS_MPI3 = $(MDEFS)-DFLASH_MPI3 $(MDEFS)-DFLASH_LIBNBC
	LIB_MPI3 = $(LIB_LIBNBC)
endif


EXE = test_pipeline

UT_FLAGS = $(MDEFS)-DUTPIPELINE_UNIT_TEST
FPPFLAGS = $(UT_FLAGS) $(FFLAGS_MPI3)

FFLAGS = $(FFLAGS_DEBUG) $(F90FLAGS) $(FPPFLAGS)
LIBS = $(LIB_MPI) $(LIB_MPI3)
LFLAGS = $(LFLAGS_DEBUG)


FSRC = UTPipeline.F90 test_pipeline.F90
FOBJ = $(FSRC:.F90=.o)


$(EXE): unlink_files link_files $(FOBJ)
	$(FCOMP) $(LFLAGS) $(EXE) $(FOBJ) $(LIBS)

link_files:
	ln -s ../UTPipeline.F90

test_pipeline.o : UTPipeline.o

$(FOBJ): %.o: %.F90
	$(FCOMP) $< $(FFLAGS) -o $@

unlink_files:
	rm -f UTPipeline.F90

clean: unlink_files
	rm -f $(FOBJ) $(EXE) *.mod core.* test*.out *__genmod.f90 valgrind.log.*

run:
	./$(EXE)
