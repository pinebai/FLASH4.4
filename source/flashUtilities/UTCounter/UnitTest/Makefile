# To build the unit test, type one of the following:
# > make dist_counter
# > make ms_counter
# dist_counter requires MPI-3; ms_counter optionally uses MPI-3
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
# > make LIBNBC=1 HEADER=Makefile.h.gnu ms_counter


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

UT_FLAGS = $(MDEFS)-DUTCOUNTER_UNIT_TEST
ms_counter: FPPFLAGS = $(UT_FLAGS) $(FFLAGS_MPI3) $(MDEFS)-DUNMET_TARGET_OK
dist_counter: FPPFLAGS = $(UT_FLAGS) $(FFLAGS_MPI3)


FFLAGS = $(FFLAGS_DEBUG) $(F90FLAGS) $(FPPFLAGS)
LIBS = $(LIB_MPI) $(LIB_MPI3)
LFLAGS = $(LFLAGS_DEBUG)
FSRC = UTCounter_sharedCounter.F90 test_shared_counter.F90
FOBJ = $(FSRC:.F90=.o)


info:
	@echo Type \'make ms_counter\' or \'make dist_counter\'

#The Master-slave counter
ms_counter: clean link_ms_counter $(FOBJ)
	$(FCOMP) $(LFLAGS) ms_counter $(FOBJ) $(LIBS)
link_ms_counter:
	@echo \*\*\*\*\*\*\*\*\*\* Linking "Master Slave" counter!!! \*\*\*\*\*\*\*\*\*\*
	ln -s ../MasterSlave/UTCounter_sharedCounter.F90

#The distributed counter
dist_counter: clean link_dist_counter $(FOBJ)
	$(FCOMP) $(LFLAGS) dist_counter $(FOBJ) $(LIBS)
link_dist_counter:
	@echo \*\*\*\*\*\*\*\*\*\* Linking "Distributed" counter!!!  DEPENDS ON MPI3!!! \*\*\*\*\*\*\*\*\*\*
	ln -s ../Distributed/UTCounter_sharedCounter.F90

UTCounter_sharedCounter.o: UTCounter_sharedCounter.F90
test_shared_counter.o: UTCounter_sharedCounter.o
$(FOBJ): %.o: %.F90
	$(FCOMP) $< $(FFLAGS) -o $@

#For clean-up
unlink_files:
	rm -f UTCounter_sharedCounter.F90
clean: unlink_files
	rm -f $(FOBJ) ms_counter dist_counter *.mod core.* test*.out *__genmod.f90 valgrind.log.*
