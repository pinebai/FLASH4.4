# To build and run the unit test, type:
# > make run
# The test should complete even though there is an invalid write to
# memory.  We can detect this invalid write using mcheck (part of
# glibc) or dmalloc.
#
# To use mcheck, type:
# > make MCHECK=1 run
# The memory violation should be detected.
#
# To use dmalloc, type:
# > make DMALLOC=1 run
# The memory violation is only detected if DMALLOC_OPTIONS is set
# appropriately.  To use good default DMALLOC_OPTIONS, type:
# > make DMALLOC=1 dmalloc_run
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

include ../../../../../sites/$(SITEDIR)/$(MAKEFILE_HEADER)

#This is the default path to the stubs.
IMPL_PATH = ..

ifdef MCHECK
	CFLAGS_EXTRA = $(MDEFS)-DMCHECK
	LIBS_EXTRA = -lmcheck
	IMPL_PATH = ../mcheck
endif


ifdef DMALLOC
	CFLAGS_EXTRA = -I$(DMALLOC_PATH)/include $(MDEFS)-DMALLOC
	LIBS_EXTRA = -L$(DMALLOC_PATH)/lib -ldmalloc
	IMPL_PATH = ../dmalloc
endif


EXE = out_of_bounds
CFLAGS = $(CFLAGS_DEBUG) $(CFLAGS_EXTRA)
LIBS = $(LIB_MPI) $(LIBS_EXTRA)
LFLAGS = $(LFLAGS_DEBUG)

SRC = out_of_bounds.c dbg_heap_check.c dbg_heap_check_impl.c
OBJ = $(SRC:.c=.o)


$(EXE): clean link_files $(OBJ)
	$(CCOMP) $(LFLAGS) $(EXE) $(OBJ) $(LIBS)

link_files:
	ln -s ../../../../flashUtilities/general/mangle_names.h
	ln -s $(IMPL_PATH)/dbg_heap_check_impl.h
	ln -s $(IMPL_PATH)/dbg_heap_check_impl.c
	ln -s ../dbg_heap_check.h
	ln -s ../dbg_heap_check.c

unlink_files:
	rm -f mangle_names.h
	rm -f dbg_heap_check.[ch]
	rm -f dbg_heap_check_impl.[ch]

out_of_bounds.o: dbg_heap_check.o
dbg_heap_check.o: dbg_heap_check_impl.o
$(OBJ): %.o: %.c
	$(CCOMP) $< $(CFLAGS) -o $@

clean: unlink_files
	rm -f $(OBJ) $(EXE) core.* dmalloc.log.*

run: $(EXE)
	./$(EXE)

dmalloc_run: $(EXE)
	env DMALLOC_OPTIONS="debug=0x4f4ed03,inter=100,log=dmalloc.log.%p" ./$(EXE)
