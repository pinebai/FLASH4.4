#	Makefile for the PARAMESH 4.0 adaptive mesh refinement library
#     [Grid/GridMain/paramesh/paramesh4/Paramesh4.0]

Grid += prolong_arrays.o  \
        mpi_amr_1blk_guardcell.o mpi_amr_1blk_guardcell_c_to_f.o \
        mpi_amr_1blk_restrict.o mpi_amr_bsort.o mpi_amr_checkpoint.o	\
        mpi_amr_comm_setup.o mpi_amr_derefine_blocks.o mpi_amr_edge_average.o \
        mpi_amr_edge_average_udt.o mpi_amr_edge_average_vdt.o	\
        mpi_amr_edge_diagonal_check.o mpi_amr_flux_conserve.o	\
        mpi_amr_flux_conserve_udt.o mpi_amr_flux_conserve_vdt.o	\
        mpi_amr_global_domain_limits.o mpi_amr_gsurr_blks.o	\
        mpi_amr_guardcell.o mpi_amr_local_surr_blks.o mpi_amr_morton.o \
        mpi_amr_morton_limits.o mpi_amr_prolong.o	\
        mpi_amr_prolong_fc_divbconsist.o mpi_amr_redist_blk.o	\
        mpi_amr_refine_blocks.o mpi_amr_refine_derefine.o mpi_amr_restrict.o \
        mpi_amr_restrict_fulltree.o \
        mpi_amr_restrict_bnd_data.o mpi_amr_restrict_bnd_data_vdt.o \
        mpi_amr_restrict_edge_data.o mpi_amr_restrict_edge_data_vdt.o \
        mpi_amr_store_comm_info.o	\
        mpi_amr_tree_setup.o mpi_get_buffer.o		         \
        mpi_get_edge_buffer.o mpi_get_flux_buffer.o mpi_lib.o	   \
        mpi_morton_bnd.o mpi_morton_bnd_fluxcon.o mpi_morton_bnd_prolong1.o \
        mpi_morton_bnd_restrict.o mpi_pack_blocks.o mpi_pack_edges.o \
        mpi_pack_fluxes.o mpi_pack_tree_info.o mpi_put_buffer.o		\
        mpi_put_edge_buffer.o mpi_put_flux_buffer.o mpi_set_message_limits.o \
        mpi_set_message_limits_unpack.o mpi_set_message_sizes.o	 \
        mpi_unpack_blocks.o mpi_unpack_edges.o mpi_unpack_fluxes.o	\
        mpi_unpack_tree_info.o rationalize_list.o mpi_amr_get_remote_block.o \
        mpi_amr_get_remote_block_fvar.o   \
        amr_1blk_bcset.o amr_1blk_cc_cp_remote.o	\
        amr_1blk_cc_prol_gen_unk_fun.o amr_1blk_cc_prol_gen_work_fun.o	\
        amr_1blk_copy_soln.o amr_1blk_ec_cp_remote.o		\
        amr_1blk_ec_prol_gen_fun.o amr_1blk_fc_cp_remote.o	\
        amr_1blk_fc_prol_gen_fun.o amr_1blk_guardcell_reset.o		\
        amr_1blk_guardcell_srl.o amr_1blk_nc_cp_remote.o	\
        amr_1blk_nc_prol_gen_fun.o amr_1blk_save_soln.o amr_1blk_t_to_perm.o \
        amr_1blk_to_perm.o amr_bcset_init.o amr_close.o	        \
        amr_initialize.o amr_mpi_find_blk_in_buffer.o amr_perm_to_1blk.o \
        amr_prolong_cc_fun_init.o amr_prolong_face_fun_init.o	\
        amr_prolong_fun_init.o amr_restrict_ec_fun.o amr_restrict_edge.o \
        amr_restrict_fc_fun.o amr_restrict_red.o amr_restrict_unk_fun.o	\
        amr_restrict_work_fun.o	\
        amr_restrict_work_fun_recip.o mesh_test.o quicksort_index.o \
        amr_system_calls.o amr_1blk_cc_prol_genorder.o	\
        amr_1blk_cc_prol_inject.o amr_1blk_cc_prol_linear.o	\
        amr_1blk_cc_prol_work_genorder.o \
        amr_1blk_cc_prol_work_mg.o \
        amr_1blk_cc_prol_work_inject.o amr_1blk_cc_prol_work_linear.o	\
        amr_1blk_fc_prol_inject.o	\
        amr_1blk_fc_prol_linear.o \
        amr_block_geometry.o \
        paramesh_mpi_interfaces.o paramesh_interfaces.o \
        mpi_wrapper_int.o mpi_wrapper_logical.o \
        mpi_wrapper_real.o amr_1blk_ec_prol_linear.o \
        amr_1blk_ec_prol_genorder.o amr_1blk_fc_prol_genorder.o	\
        amr_1blk_nc_prol_linear.o amr_1blk_nc_prol_genorder.o     \
        mpi_amr_boundary_block_info.o mpi_mort_comm_for_surrblks.o  \
        amr_set_runtime_parameters.o      \
        amr_1blk_fc_clean_divb.o physicaldata.o \
        mpi_amr_checkpoint_default.o \
        mpi_amr_checkpoint_hdf5.o \
        mpi_amr_checkpoint_mpiio.o \
        mpi_amr_plotfile_chombo.o \
        mpi_amr_singular_line.o \
        mpi_amr_test_neigh_values.o \
        mpi_put_edge_buffer_1blk.o \
        amr_1blk_cc_prol_user.o \
        amr_1blk_cc_prol_work_user.o \
        amr_1blk_ec_prol_user.o \
        amr_1blk_fc_prol_user.o \
        amr_1blk_nc_prol_user.o \
        amr_restrict_ec_genorder.o \
        amr_restrict_ec_user.o \
        amr_restrict_fc_genorder.o \
        amr_restrict_fc_ins.o \
        amr_restrict_fc_user.o \
        amr_restrict_nc_fun.o \
        amr_restrict_nc_genorder.o \
        amr_restrict_nc_user.o \
        amr_restrict_unk_genorder.o \
        amr_restrict_unk_user.o \
        amr_restrict_work_genorder.o \
        amr_restrict_work_user.o \
        amr_reorder_grid.o \
        user_coord_transfm.o \
        mpi_wrapper_dble.o\
	amr_1blk_fc_prol_dbz.o \
	paramesh_comm_data.o


amr_1blk_cc_cp_remote.o \
	amr_1blk_cc_prol_genorder.o \
	amr_1blk_cc_prol_work_genorder.o \
	amr_1blk_cc_prol_work_mg.o \
	amr_1blk_copy_soln.o \
	amr_1blk_fc_cp_remote.o \
	amr_close.o \
	amr_mpi_find_blk_in_buffer.o \
	amr_reorder_grid.o \
	amr_set_runtime_parameters.o: Flash.h paramesh_reorder.h clean_divb.fh
mpi_amr_1blk_guardcell_c_to_f.o \
	mpi_amr_1blk_guardcell.o \
	mpi_amr_1blk_restrict.o \
	mpi_amr_checkpoint_default.o \
	mpi_amr_comm_setup.o \
	mpi_amr_derefine_blocks.o \
	mpi_amr_edge_average_udt.o \
	mpi_amr_edge_diagonal_check.o \
	mpi_amr_get_remote_block.o \
	mpi_amr_get_remote_block_fvar.o \
	mpi_amr_gsurr_blks.o \
	mpi_amr_guardcell.o \
	mpi_amr_local_surr_blks.o \
	mpi_amr_mirror_blks.o \
	mpi_amr_mirror_solution.o \
	mpi_amr_morton.o \
	mpi_amr_morton_limits.o \
	mpi_amr_prolong.o \
	mpi_amr_prolong_fc_divbconsist.o \
	mpi_amr_redist_blk.o \
	mpi_amr_refine_blocks.o \
	mpi_amr_refine_derefine.o \
	mpi_amr_restrict_bnd_data.o \
	mpi_amr_restrict_edge_data.o \
	mpi_amr_restrict_edge_data_vdt.o \
	mpi_amr_restrict.o \
	mpi_amr_restrict_fulltree.o \
	mpi_amr_singular_line.o \
	mpi_amr_test_neigh_values.o \
	mpi_amr_tree_setup.o \
	mpi_get_buffer.o \
	mpi_get_edge_buffer.o \
	mpi_get_flux_buffer.o \
	mpi_lib.o \
	mpi_mort_comm_for_surrblks.o \
	mpi_morton_bnd.o \
	mpi_morton_bnd_fluxcon.o \
	mpi_morton_bnd_prolong1.o \
	mpi_morton_bnd_restrict.o \
	mpi_morton.o \
	mpi_pack_blocks.o \
	mpi_pack_edges.o \
	mpi_pack_fluxes.o \
	mpi_pack_tree_info.o \
	mpi_put_buffer.o \
	mpi_put_edge_buffer.o \
	mpi_put_flux_buffer.o \
	mpi_set_message_sizes.o \
	mpi_unpack_blocks.o \
	mpi_unpack_edges.o \
	mpi_unpack_fluxes.o \
	mpi_unpack_tree_info.o: Flash.h paramesh_reorder.h clean_divb.fh
amr_initialize.o gr_amr_dump_runtime_parameters.o \
	io.o \
	paramesh_dimensions.o \
	physicaldata.o \
	prolong_arrays.o \
	timings.o \
	tree.o \
	workspace.o : Flash.h paramesh_reorder.h clean_divb.fh
#mesh_test.o paramesh_interfaces.o parse_solver.o: Flash.h paramesh_reorder.h clean_divb.fh
ifdef MODUPPERCASE
IO.mod \
	PARAMESH_DIMENSIONS.mod \
	PHYSICALDATA.mod \
	PROLONG_ARRAYS.mod \
	TIMINGS.mod \
	TREE.mod \
	MPI_MORTON.mod \
	WORKSPACE.mod : Flash.h paramesh_reorder.h clean_divb.fh
else
io.mod \
	paramesh_dimensions.mod \
	physicaldata.mod \
	prolong_arrays.mod \
	timings.mod \
	tree.mod \
	mpi_morton.mod \
	workspace.mod : Flash.h paramesh_reorder.h clean_divb.fh
endif

mpi_amr_checkpoint.o : io.o

Grid += tree.o paramesh_dimensions.o \
        workspace.o mpi_morton.o timings.o io.o constants.o

Grid += gr_amr_dump_runtime_parameters.o

Grid += flash_convert_cc_hook.o flash_unconvert_cc_hook.o \
	gr_sanitizeDataAfterInterp.o gr_flashHook_interfaces.o gr_commSetUp.o \
	gr_packBCs.o \
	gr_setGcFillNLayers.o gr_makeMaskConsistent.o gr_setMasks.o \
	gr_findNeghID.o gr_findAllNeghID.o gr_checkGridState.o

Grid += amr_get_new_nodetypes.o amr_mg_init.o amr_mg_prolong.o \
	mpi_amr_store_comm_info_mg.o amr_mg_common.o \
	amr_mg_morton_process.o amr_mg_restrict.o

