#ifndef CHOMBO_F_C_API_H
#define CHOMBO_F_C_API_H

#include <cstdlib>
#include <cassert>
#include "constants.h"
#include "Flash.h"
#include "flash_ctypes.h"

// The functions in this file are callable by Fortran and are here
// to operate on our Chombo grid object
#ifdef FLASH_GRID_UG
#define CHOMBO_UG
#include "chombo_uniform_grid.h"
Chombo_Uniform_Grid *pChomboGrid = 0;
#else
#define CHOMBO_AMR
#include "chombo_adaptive_grid.h"
Chombo_Adaptive_Grid *pChomboGrid = 0;
#endif

extern "C" void c_free(void **p);

// We should really have a header file for Driver_abortFlashC.  I did not
// add one because it means touching a huge number of files.
extern "C" int Driver_abortFlashC(char* message);

extern "C" void ch_define_uniform_grid(const flash_ug_info_t *pFlashUGInfo,
				       const int meshStringLens[],
				       const char meshStrings[],
                       bool restart);

extern "C" void ch_define_adaptive_grid(const flash_amr_info_t *pFlashAMRInfo,
					const int meshStringLens[],
					const char meshStrings[],
                    bool restart);

extern "C" void ch_is_initial_refinement_done(int *pRefComplete);

extern "C" void ch_build_initial_grid();

extern "C" void ch_refine_initial_grid();

extern "C" void ch_finalize_initial_grid();

extern "C" void ch_regrid(const int baseLevel);

extern "C" void ch_restrict_all_levels();

extern "C" void ch_fill_guardcells();

extern "C" void ch_write_checkpoint(
    const char filename[],
	double simTime,
	double dt,
    const named_vals_t *scalars,
    const named_vals_t *runparms);

extern "C" void ch_read_checkpoint(
    const char filename[],
	double *simTime,
	double *dt,
    named_vals_t *scalars,
    named_vals_t *runparms);

extern "C" void ch_get_blk_ptr(const int blkID,
			       const int gridDataStruct,
			       void **ppv);

extern "C" void ch_get_block_ids(const int blkType,
				 const int level,
				 int blkIDs[],
				 int *pNumBlks);

extern "C" void ch_get_num_blocks(int *pNumBlks);

extern "C" void ch_get_cell_coords(const int blkID,
				   const int axis,
				   const int edge,
				   const int size,
				   const int guardcell,
				   double coords[]);

extern "C" void ch_get_single_cell_coords(const int blkID,
					  const int edge,
					  const int beginCount,
					  int ind[],
					  double coords[]);

extern "C" void ch_get_box_info(const int blkID,
				const int gridDataStruct,
				box_info_t *pBoxInfo);

extern "C" void ch_zero_flux_data();

extern "C" void ch_put_flux_data(const void *pFluxes,
				 const double dt,
				 const int blkID,
				 const int axis);

extern "C" void ch_reflux();

extern "C" void ch_post_time_step();

extern "C" void ch_finalize();

// Helper functions below this point (i.e. not callable by Fortran).

void assert_valid_simulation();

void extract_mesh_strings(const int meshNumVars[],
			  const int meshTypes[],
			  const int meshStringLens[],
			  const char meshStrings[],
			  mesh_info_t &meshInfo);

void unit_base_box_info(box_info_t *pBoxInfo);

#endif
