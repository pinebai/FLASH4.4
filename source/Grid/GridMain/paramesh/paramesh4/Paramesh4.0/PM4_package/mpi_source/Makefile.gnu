sources := \
mpi_amr_singular_line.F90 \
mpi_amr_block_boundary_tecplot.F90 \
mpi_amr_1blk_guardcell.F90 \
mpi_amr_1blk_guardcell_c_to_f.F90 \
mpi_amr_1blk_restrict.F90 \
mpi_amr_bsort.F90 \
mpi_amr_comm_setup.F90 \
mpi_amr_edge_average.F90 \
mpi_amr_edge_average_udt.F90 \
mpi_amr_edge_average_vdt.F90 \
mpi_amr_edge_diagonal_check.F90 \
mpi_amr_flux_conserve.F90 \
mpi_amr_flux_conserve_udt.F90 \
mpi_amr_flux_conserve_vdt.F90 \
mpi_amr_get_remote_block.F90 \
mpi_amr_get_remote_block_fvar.F90 \
mpi_amr_global_domain_limits.F90 \
mpi_amr_gsurr_blks.F90 \
mpi_amr_guardcell.F90 \
mpi_amr_local_surr_blks.F90 \
mpi_amr_morton_limits.F90 \
mpi_amr_prolong.F90 \
mpi_amr_prolong_fc_divbconsist.F90 \
mpi_amr_refine_derefine.F90 \
mpi_amr_restrict.F90 \
mpi_amr_restrict_fulltree.F90 \
mpi_amr_restrict_bnd_data_vdt.F90 \
mpi_amr_restrict_edge_data_vdt.F90 \
mpi_amr_store_comm_info.F90 \
mpi_amr_timing_report.F90 \
mpi_amr_tree_setup.F90 \
mpi_get_buffer.F90 \
mpi_get_edge_buffer.F90 \
mpi_get_flux_buffer.F90 \
mpi_amr_mirror_blks.F90 \
mpi_amr_mirror_solution.F90 \
mpi_mort_comm_for_surrblks.F90 \
mpi_morton_bnd.F90 \
mpi_morton_bnd_fluxcon.F90 \
mpi_morton_bnd_prolong1.F90 \
mpi_morton_bnd_restrict.F90 \
mpi_pack_blocks.F90 \
mpi_unpack_blocks.F90 \
mpi_put_buffer.F90 \
mpi_pack_edges.F90 \
mpi_pack_fluxes.F90 \
mpi_put_edge_buffer.F90 \
mpi_put_edge_buffer_1blk.F90 \
mpi_put_flux_buffer.F90 \
mpi_set_message_limits.F90 \
mpi_set_message_sizes.F90 \
mpi_unpack_edges.F90 \
mpi_unpack_fluxes.F90 \
rationalize_list.F90 \
mpi_amr_checkpoint.F90 \
mpi_amr_checkpoint_default.F90 \
mpi_amr_checkpoint_hdf5.F90 \
mpi_amr_checkpoint_mpiio.F90 \
read_blocks_hdf5_r4.c \
read_blocks_hdf5_r8.c \
write_blocks_hdf5_r4.c \
write_blocks_hdf5_r8.c \
mpi_amr_plotfile_chombo.F90 \
write_blocks_chombo_r4.c \
write_blocks_chombo_r8.c \
mpi_amr_derefine_blocks.F90 \
mpi_amr_morton.F90 \
mpi_amr_redist_blk.F90 \
mpi_amr_refine_blocks.F90 \
mpi_amr_restrict_bnd_data.F90 \
mpi_amr_restrict_edge_data.F90 \
mpi_amr_boundary_block_info.F90 \
mpi_amr_get_bc_settings.F90 \
mpi_amr_test_neigh_values.F90 \
mpi_lib.F90 \
mpi_pack_tree_info.F90 \
mpi_unpack_tree_info.F90 \
mpi_wrapper_int.F90 \
mpi_wrapper_logical.F90 \
mpi_wrapper_real.F90 \
mpi_wrapper_dble.F90

%.o:%.F90
	$(FC) -c $(FFLAGS) $<

objects := $(sources:.F90=.o)
objects := $(objects:.c=.o)

vpath %.fh ../headers

libmpi_paramesh.a: $(objects)
	$(AR) $(ARFLAGS) $@ $^
	cp -f libmpi_paramesh.a ../libs/.

ifdef MY_CPP
GNUmakefile.include: $(sources)
	find . -name \*.F90 | xargs $(MY_CPP) > $@
include GNUmakefile.include
else
$(objects): $(wildcard *.fh)
endif

.PHONY: clean
clean:
	$(RM) libmpi_paramesh.a *.o *.i *~ GNUmakefile.include
