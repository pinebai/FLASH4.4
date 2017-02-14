# This Makefile will use the macro definitions from the
# FLASH Makefile.h in the sites directory.  The site
# directory name is obtained from the hostname command.
#
# To build the unit test:
# > make
#
# If the sites directory has a different name to the return
# value of hostname command then you must define SITE macro to be
# the sites directory of interest:
# > make SITE=code.uchicago.edu
#

# Include the appropriate Makefile.h.
ifdef SITE
	SITEDIR=$(SITE)
else
	SITEDIR:=$(shell hostname -f)
endif
include ../../../../../sites/$(SITEDIR)/Makefile.h


EXE = ut_sysMemUnitTest1

FFLAGS = $(FFLAGS_OPT) $(F90FLAGS)
CFLAGS = $(CFLAGS_OPT) $(CFLAGS_MPI)
LIBS = $(LIB_MPI)
LFLAGS = $(LFLAGS_OPT)

SUPPORT_HEADERS = constants.h \
	ut_sysMemCTypes.h ut_sysMemMallinfo.h ut_sysMemProc.h \
	ut_sysMemBGKernel.h ut_sysMemRusage.h ut_sysMemAutoDetect.h ut_sysMem.h
SUPPORT_CFILES = ut_sysMemMallinfo.c ut_sysMemProc.c ut_sysMemBGKernel.c \
	ut_sysMemRusage.c ut_sysMemAutoDetect.c
SUPPORT_FFILES = ut_sysMemData.F90 ut_sysMemInterface.F90 \
	ut_sysMemCData.F90 ut_sysMemCInterface.F90 ut_sysMemStats.F90

CSRC = $(SUPPORT_CFILES)
COBJ = $(CSRC:.c=.o)

FSRC = $(SUPPORT_FFILES) ut_sysMemUnitTest1.F90
FOBJ = $(FSRC:.F90=.o)


$(EXE): unlink_files link_files $(COBJ) $(FOBJ)
	$(FCOMP) $(LFLAGS) $(EXE) $(COBJ) $(FOBJ) $(LIBS)

link_files:
	ln -s ../../../../Simulation/constants.h
	ln -s ../ut_sysMem.h
	ln -s ../ut_sysMemData.F90
	ln -s ../ut_sysMemInterface.F90
	ln -s ../f2003/ut_sysMemCData.F90
	ln -s ../f2003/ut_sysMemCInterface.F90
	ln -s ../f2003/ut_sysMemAutoDetect.c
	ln -s ../f2003/ut_sysMemAutoDetect.h
	ln -s ../f2003/ut_sysMemCTypes.h
	ln -s ../f2003/ut_sysMemBGKernel.c
	ln -s ../f2003/ut_sysMemBGKernel.h
	ln -s ../f2003/ut_sysMemMallinfo.c
	ln -s ../f2003/ut_sysMemMallinfo.h
	ln -s ../f2003/ut_sysMemProc.c
	ln -s ../f2003/ut_sysMemProc.h
	ln -s ../f2003/ut_sysMemRusage.c
	ln -s ../f2003/ut_sysMemRusage.h
	ln -s ../f2003/ut_sysMemStats.F90

$(COBJ): %.o: %.c
	$(CCOMP) $< $(CFLAGS) $(CPPFLAGS) -o $@
$(FOBJ): %.o: %.F90
	$(FCOMP) $< $(FFLAGS) $(FPPFLAGS) -o $@

unlink_files:
	rm -f $(SUPPORT_HEADERS) $(SUPPORT_CFILES) $(SUPPORT_FFILES)

clean: unlink_files
	rm -f $(COBJ) $(FOBJ) $(EXE) *.mod

run:
	./$(EXE)
