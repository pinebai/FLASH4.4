# This Makefile will use the macro definitions from the
# FLASH Makefile.h.chombo in the sites directory.  The site
# directory name is obtained from the hostname command.
#
# To build the unit test:
# > make all
#
# If the sites directory has a different name to the return
# value of hostname command then you must define SITE macro to be
# the sites directory of interest:
# > make all SITE=code.uchicago.edu
#

# Include the appropriate Makefile.h.
ifdef SITE
	SITEDIR=$(SITE)
else
	SITEDIR:=$(shell hostname)
endif
include ../../../../../../../sites/$(SITEDIR)/Makefile.h.chombo


EXE = unit_test_2

# These macros are obtained from your Makefile.h.chombo.
#   We do this because CFLAGS_CHOMBO and LIB_CHOMBO depend
#   heavily on how you build your Chombo library
#   (see sites/code.uchicago.edu/Makefile.h.chombo).
FFLAGS = $(FFLAGS_DEBUG)
CFLAGS = $(CFLAGS_DEBUG) $(CFLAGS_CHOMBO) $(CFLAGS_HDF5) $(CFLAGS_MPI)
LFLAGS = $(LFLAGS_DEBUG)
LIB = $(LIB_CHOMBO) $(LIB_HDF5) $(LIB_MPI)


CSRC = Driver_abortFlashC.c
COBJ = $(CSRC:.c=.o)
CXXSRC = chombo_f_c_api.C chombo_uniform_grid.C
CXXOBJ = $(CXXSRC:.C=.o)
FSRC     = Driver_interface.F90 flash_interfaces.F90 \
	flash_ftypes.F90 chombo_f_c_interface.F90 \
	Grid_getBlkIndexLimits.F90 Grid_getBlkPtr.F90 \
	flash_subroutines.F90 main.F90
FOBJ = $(FSRC:.F90=.o)


$(CXXOBJ): %.o: %.C
	$(CPPCOMP) -c $< $(CFLAGS) -o $@

$(COBJ): %.o: %.c
	$(CCOMP) -c $< $(CFLAGS) -o $@

$(FOBJ): %.o: %.F90
	$(FCOMP) -c $< $(FFLAGS) -o $@


chombo_f_c_api.o: chombo_f_c_api.C chombo_f_c_api.h \
	chombo_uniform_grid.h \
	Flash.h constants.h flash_ctypes.h
chombo_uniform_grid.o: chombo_uniform_grid.C chombo_uniform_grid.h \
	flash_ctypes.h constants.h flash_bool.h

flash_subroutines.o: chombo_f_c_interface.F90 flash_interfaces.F90 \
	Driver_interface.F90 Flash.h constants.h
chombo_f_c_interface.o: chombo_f_c_interface.F90 flash_ftypes.F90 \
	Flash.h constants.h

flash_interfaces.o : flash_interfaces.F90 Flash.h constants.h
Grid_getBlkPtr.o : Grid_getBlkPtr.F90 flash_ftypes.F90 \
	Driver_interface.F90 Flash.h constants.h
Grid_getBlkIndexLimits.o : Grid_getBlkIndexLimits.F90 flash_ftypes.F90 \
	Driver_interface.F90 Flash.h constants.h


all: $(EXE) run

$(EXE): link_files $(CXXOBJ) $(COBJ) $(FOBJ)
	$(FCOMP) $(LFLAGS) $(EXE) $(CXXOBJ) $(COBJ) $(FOBJ) $(LIB)

link_files:
	ln -s ../../../Grid_getBlkPtr.F90
	ln -s ../../../Grid_getBlkIndexLimits.F90
	ln -s ../../../../../../Simulation/constants.h
	ln -s ../../../../../../flashUtilities/general/mangle_names.h
	ln -s ../../../../../../Driver/DriverMain/Driver_abortFlashC.c
	ln -s ../../chombo_f_c_api.C
	ln -s ../../chombo_f_c_api.h
	ln -s ../../chombo_f_c_interface.F90
	ln -s ../../chombo_uniform_grid.C
	ln -s ../../chombo_uniform_grid.h
	ln -s ../../flash_ctypes.h
	ln -s ../../flash_bool.h
	ln -s ../../flash_ftypes.F90

run:
	mpirun -np 4 ./$(EXE)

unlink_files:
	rm -f Grid_getBlkPtr.F90
	rm -f Grid_getBlkIndexLimits.F90
	rm -f constants.h
	rm -f mangle_names.h
	rm -f Driver_abortFlashC.c
	rm -f chombo_f_c_api.C
	rm -f chombo_f_c_api.h
	rm -f chombo_f_c_interface.F90
	rm -f chombo_uniform_grid.C
	rm -f chombo_uniform_grid.h
	rm -f flash_ctypes.h
	rm -f flash_bool.h
	rm -f flash_ftypes.F90

clean: unlink_files
	rm -f $(CXXOBJ) $(COBJ) $(FOBJ) $(EXE) *.mod *__genmod.* pout.[0-3]
