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
include ../../../../../../../sites/$(SITEDIR)/Makefile.h


EXE = Test_parallel_write

CPPFLAGS =
CFLAGS = $(CFLAGS_OPT) $(CFLAGS_NCMPI) $(CFLAGS_MPI)
LIBS = $(LIB_NCMPI) $(LIB_MPI)
LFLAGS = $(LFLAGS_OPT)

SUPPORT_HEADERS = constants.h io_flash.h \
	io_ncmpi_type.h io_ncmpi_attribute.h io_use_grid_mpi_types.h \
	io_set_grid_mpi_types.h io_ncmpi_xfer_mesh_dataset.h \
	io_ncmpi_nonblocking.h io_ncmpi_create_dataset.h \
	mangle_names.h
SUPPORT_FILES = Driver_abortFlashC.c io_ncmpi_type.c io_ncmpi_attribute.c \
	io_use_grid_mpi_types.c io_set_grid_mpi_types.c \
	io_ncmpi_xfer_mesh_dataset.c io_ncmpi_nonblocking.c \
	io_ncmpi_create_dataset.c
CSRC = $(SUPPORT_FILES) Test_parallel_write.c
COBJ = $(CSRC:.c=.o)


$(EXE): link_files $(COBJ)
	$(CCOMP) $(LFLAGS) $(EXE) $(COBJ) $(LIBS)

link_files:
	ln -s ../../../../hdf5/parallel/PM_argonne/UnitTest/Driver_abortFlashC.c
	ln -s ../../../../hdf5/parallel/PM_argonne/UnitTest/constants.h
	ln -s ../../../../io_flash.h
	ln -s ../../../io_ncmpi_type.c
	ln -s ../../../io_ncmpi_type.h
	ln -s ../../../io_ncmpi_attribute.c
	ln -s ../../../io_ncmpi_attribute.h
	ln -s ../../../../../IOTypes/io_use_grid_mpi_types.c
	ln -s ../../../../../IOTypes/io_use_grid_mpi_types.h
	ln -s ../../../../../IOTypes/io_set_grid_mpi_types.c
	ln -s ../../../../../IOTypes/io_set_grid_mpi_types.h
	ln -s ../../../../../IOTypes/pnetcdf/io_ncmpi_xfer_mesh_dataset.c	
	ln -s ../../../../../IOTypes/pnetcdf/io_ncmpi_xfer_mesh_dataset.h
	ln -s ../../../../../IOTypes/pnetcdf/io_ncmpi_nonblocking.c
	ln -s ../../../../../IOTypes/pnetcdf/io_ncmpi_nonblocking.h
	ln -s ../../../io_ncmpi_create_dataset.c
	ln -s ../../../io_ncmpi_create_dataset.h
	ln -s ../../../../../../flashUtilities/general/mangle_names.h

$(COBJ): %.o: %.c
	$(CCOMP) $< $(CFLAGS) $(CPPFLAGS) -o $@

unlink_files:
	rm -f $(SUPPORT_HEADERS) $(SUPPORT_FILES)

clean: unlink_files
	rm -f $(COBJ) $(EXE) \
	Test_parallel_write.pnetcdf pnetcdf_new.out

run:
	mpirun -np 4 ./$(EXE)
	ncmpidump Test_parallel_write.pnetcdf > pnetcdf_new.out
