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


EXE = Test_parallel_write

#These pre-processor flags are required in some SUPPORT_FILES.
CPPFLAGS = -DIO_HDF5_PARALLEL -DN_DIM=1

#For IBM
#CPPFLAGS = -DIBM -DIO_HDF5_PARALLEL -DN_DIM=1

CFLAGS = $(CFLAGS_OPT) $(CFLAGS_HDF5) $(CFLAGS_MPI)
LIBS = $(LIB_HDF5) $(LIB_MPI)
LFLAGS = $(LFLAGS_OPT)

SUPPORT_HEADERS = Flash.h constants.h io_flash.h hdf5_flash.h \
	io_h5create_dataset.h io_h5_attribute.h io_h5_type.h \
	io_h5_xfer_wrapper.h io_h5_type_matched_xfer.h io_h5_xfer.h
SUPPORT_FILES = Driver_abortFlashC.c io_h5_type.c io_h5create_dataset.c \
	io_h5_attribute.c io_h5_type_matched_xfer.c \
	io_h5_xfer_wrapper.c io_h5_xfer.c
CSRC = $(SUPPORT_FILES) Test_parallel_write.c
COBJ = $(CSRC:.c=.o)


$(EXE): link_files $(COBJ)
	$(CCOMP) $(LFLAGS) $(EXE) $(COBJ) $(LIBS)

link_files:
	ln -s ../Driver_abortFlashC.c
	ln -s ../Flash.h
	ln -s ../constants.h
	ln -s ../../../../../io_flash.h
	ln -s ../../../../hdf5_flash.h
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

$(COBJ): %.o: %.c
	$(CCOMP) $< $(CFLAGS) $(CPPFLAGS) -o $@

unlink_files:
	rm -f $(SUPPORT_HEADERS) $(SUPPORT_FILES)

clean: unlink_files
	rm -f $(COBJ) $(EXE) \
	CheckpointFile.hdf5 PlotFile.hdf5 ch5_new.out ph5_new.out

run:
	mpirun -np 4 ./$(EXE)
	h5dump CheckpointFile.hdf5 > ch5_new.out
	h5dump PlotFile.hdf5 > ph5_new.out
