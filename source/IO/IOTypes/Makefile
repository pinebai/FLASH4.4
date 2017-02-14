IO += io_typeData.o io_c_type_interface.o io_set_grid_mpi_types.o io_use_grid_mpi_types.o io_xfer_mesh_data.o io_xfer_mesh_dataset.o io_create_grid_header.o io_getZeroBasedBlkSubarray.o io_getZeroBasedVarInfo.o io_typeInit.o io_typeFinalize.o io_do_xfer.o

io_xfer_mesh_dataset.o : constants.h mangle_names.h io_flash.h Flash.h

# The following GNU make special prevents that apparent dependencies
# on the file iso_c_binding.mod, which does usually not actually exist
# in the object directory but refers to the ISO_C_BINDING module known
# to FORTRAN compilers internally, trigger unnecessary recompilation
# of files that refer to the ISO_C_BINDING module.
.SECONDARY: iso_c_binding.mod
