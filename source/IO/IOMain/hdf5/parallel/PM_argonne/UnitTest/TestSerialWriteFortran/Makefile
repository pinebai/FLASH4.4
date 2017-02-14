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
	SITEDIR:=$(shell hostname)
endif
include ../../../../../../../../sites/$(SITEDIR)/Makefile.h


EXE = Test_write_fortran

#These pre-processor flags are required in some SUPPORT_FILES.
CPPFLAGS = -DIO_HDF5_PARALLEL -DFLASH_IO_HDF5 -DN_DIM=1
FPPFLAGS = $(CPPFLAGS)

#For IBM
#CPPFLAGS = -DIBM -DIO_HDF5_PARALLEL -DFLASH_IO_HDF5 -DN_DIM=1
#FPPFLAGS = -Wl,-DIBM -Wl,-DFLASH_IO_HDF5 -Wl,-DH5_USE_16_API

FFLAGS = $(FFLAGS_OPT) $(F90FLAGS)
CFLAGS = $(CFLAGS_OPT) $(CFLAGS_HDF5) $(CFLAGS_MPI)
LIBS = $(LIB_HDF5) $(LIB_MPI)
LFLAGS = $(LFLAGS_OPT)

SUPPORT_HEADERS = Flash.h constants.h io_flash.h hdf5_flash.h \
	io_h5create_dataset.h io_h5_attribute.h io_h5_type.h \
	io_h5_xfer_wrapper.h io_h5_type_matched_xfer.h io_h5_xfer.h \
	io_attribute.h io_create_dataset.h io_xfer_cont_slab.h \
	mangle_names.h
SUPPORT_FILES = Driver_abortFlashC.c io_h5_type.c io_h5create_dataset.c \
	io_h5_attribute.c io_h5_type_matched_xfer.c \
	io_h5_xfer_wrapper.c io_h5_xfer.c \
	io_attribute.c io_create_dataset.c \
	io_xfer_cont_slab.c
CSRC = $(SUPPORT_FILES) set_hdf5_c.c
COBJ = $(CSRC:.c=.o)

FSRC = Test_write_fortran.F90
FOBJ = $(FSRC:.F90=.o)


$(EXE): link_files $(COBJ) $(FOBJ)
	$(FCOMP) $(LFLAGS) $(EXE) $(COBJ) $(FOBJ) $(LIBS)

link_files:
	ln -s ../Driver_abortFlashC.c
	ln -s ../Flash.h
	ln -s ../constants.h
	ln -s ../../../../../io_flash.h
	ln -s ../../../../hdf5_flash.h
	ln -s ../../../../../../../flashUtilities/general/mangle_names.h
	ln -s ../../../../io_h5create_dataset.c
	ln -s ../../../../io_h5create_dataset.h
	ln -s ../../../../io_h5_attribute.c
	ln -s ../../../../io_h5_attribute.h
	ln -s ../../../../io_h5_type.c
	ln -s ../../../../io_h5_type.h
	ln -s ../../../../io_h5_xfer_wrapper.c
	ln -s ../../../../io_h5_xfer_wrapper.h
	ln -s ../../../../io_h5_type_matched_xfer.c
	ln -s ../../../../io_h5_type_matched_xfer.h
	ln -s ../../../../io_h5_xfer.c
	ln -s ../../../../io_h5_xfer.h
	ln -s ../../../../../io_attribute.c
	ln -s ../../../../../io_attribute.h
	ln -s ../../../../../io_create_dataset.c
	ln -s ../../../../../io_create_dataset.h
	ln -s ../../../../../io_xfer_cont_slab.c
	ln -s ../../../../../io_xfer_cont_slab.h

$(COBJ): %.o: %.c
	$(CCOMP) $< $(CFLAGS) $(CPPFLAGS) -o $@
$(FOBJ): %.o: %.F90
	$(FCOMP) $< $(FFLAGS) $(FPPFLAGS) -o $@

unlink_files:
	rm -f $(SUPPORT_HEADERS) $(SUPPORT_FILES)

clean: unlink_files
	rm -f $(COBJ) $(FOBJ) $(EXE) Test_write_fortran.hdf5 h5_new.out

run:
	./$(EXE)
	h5dump Test_write_fortran.hdf5 > h5_new.out
