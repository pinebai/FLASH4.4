IO += IO_data.o \
	io_getAllScalars.o io_getNumScalars.o \
	io_prepareListsWrite.o io_prepareListsRead.o \
	io_finalizeListsRead.o io_finalizeListsWrite.o \
	io_bcastScalars.o io_getOutputName.o \
	io_prepareSimInfo.o \
	io_isPlotVar.o \
	io_setPrevScalar.o \
	io_getVarExtrema.o \
	io_memoryReport.o \
	io_restrictBeforeWrite.o \
	io_mpi_type.o \
	io_xfer_cont_slab.o \
	io_xfer_tree_data.o \
	io_attribute.o \
	io_create_dataset.o \
	io_c_interface.o \
	io_repack_data.o \
	io_checkBlockShape.o

io_flash.h : Flash.h; test -f $@ && ${TOUCH} $@

IO_output.o : IO_data.o Grid_sendOutputData.o
IO_init.o : RuntimeParameters_get.o IO_data.o
IO_data.o : nameValueLL_data.o
IO_setScalar.o : IO_data.o
IO_getScalar.o : IO_data.o
io_getNumScalars.o : IO_data.o
io_getAllScalars.o : IO_data.o
io_getVarExtrema.o : IO_data.o
io_bcastScalars.o : IO_data.o
IO_writePlotfile.o : io_getOutputName.o
IO_writeCheckpoint.o : Logfile_stamp.o
