!!****ih* source/Grid/GridMain/Chombo/wrapper/chombo_f_c_interface
!!
!! NAME
!!
!!  chombo_f_c_interface
!!
!! SYNOPSIS
!!
!!  use chombo_f_c_interface
!!
!! DESCRIPTION
!!
!! These are the interfaces for our C++ functions that interact
!! with Chombo library.
!! 
!!***

#include "constants.h"
#include "Flash.h"

module chombo_f_c_interface
  implicit none

  interface
     subroutine c_free(p) bind(c,name='c_free')
       use iso_c_binding
       type(c_ptr), intent(inout) :: p
     end subroutine
  end interface
  
  interface
     subroutine ch_define_adaptive_grid &
          (flashAMRInfo,meshStringLens,meshStrings,restart) &
          bind(c,name='ch_define_adaptive_grid')
       use iso_c_binding
       use flash_ftypes, only : flash_amr_info_t
       type(flash_amr_info_t), intent(in) :: flashAMRInfo
       integer(c_int), dimension(*), intent(in) :: meshStringLens
       character(kind=c_char), dimension(*), intent(in) :: meshStrings
       logical(c_bool), value, intent(in) :: restart
     end subroutine ch_define_adaptive_grid
  end interface

  interface
     subroutine ch_define_uniform_grid &
          (flashUGInfo,meshStringLens,meshStrings,restart) &
          bind(c,name='ch_define_uniform_grid')
       use iso_c_binding
       use flash_ftypes, only : flash_ug_info_t
       type(flash_ug_info_t), intent(in) :: flashUGInfo
       integer(c_int), dimension(*), intent(in) :: meshStringLens
       character(kind=c_char), dimension(*), intent(in) :: meshStrings
       logical(c_bool), value, intent(in) :: restart
     end subroutine ch_define_uniform_grid
  end interface

  interface
     subroutine ch_is_initial_refinement_done &
          (pIsRefComplete) &
          bind(c,name="ch_is_initial_refinement_done")
       use iso_c_binding, only : c_int
       integer(c_int), intent(out) :: pIsRefComplete
     end subroutine ch_is_initial_refinement_done
  end interface
  
  interface
     subroutine ch_build_initial_grid &
          () &
          bind(c,name='ch_build_initial_grid')
     end subroutine ch_build_initial_grid
  end interface

  interface
     subroutine ch_refine_initial_grid &
          () &
          bind(c,name='ch_refine_initial_grid')
     end subroutine ch_refine_initial_grid
  end interface

  interface
     subroutine ch_finalize_initial_grid &
          () &
          bind(c,name='ch_finalize_initial_grid')
     end subroutine ch_finalize_initial_grid
  end interface

  interface
     subroutine ch_finalize &
          () &
          bind(c,name='ch_finalize')
     end subroutine ch_finalize
  end interface

  interface
     subroutine ch_regrid &
          (baseLevel) &
          bind(c,name='ch_regrid')
       use iso_c_binding, only : c_int
       integer(c_int), value, intent(in) :: baseLevel
     end subroutine ch_regrid
  end interface

  interface
     subroutine ch_restrict_all_levels &
          () &
          bind(c,name='ch_restrict_all_levels')
     end subroutine ch_restrict_all_levels
  end interface

  interface
     subroutine ch_fill_guardcells &
          () &
          bind(c,name='ch_fill_guardcells')
     end subroutine ch_fill_guardcells
  end interface

  interface
     subroutine ch_write_checkpoint(filename,simTime,dt,scalars,runparms) &
          bind(c,name='ch_write_checkpoint')
       use iso_c_binding!, only : c_char, c_double, c_ptr
       use flash_ftypes, only: named_vals_t
       character(len=1,kind=c_char), dimension(*), intent(in) :: filename
       real(c_double), value, intent(in) :: simTime
       real(c_double), value, intent(in) :: dt
       type(named_vals_t), intent(in) :: scalars, runparms
     end subroutine ch_write_checkpoint
  end interface
  
  interface
     subroutine ch_read_checkpoint(filename,simTime,dt,scalars,runparms) &
          bind(c,name='ch_read_checkpoint')
       use iso_c_binding!, only : c_char, c_double, c_ptr
       use flash_ftypes, only: named_vals_t
       character(len=1,kind=c_char), dimension(*), intent(in) :: filename
       real(c_double), intent(out) :: simTime
       real(c_double), intent(out) :: dt
       type(named_vals_t), intent(out) :: scalars, runparms
     end subroutine ch_read_checkpoint
  end interface
  
  interface
     subroutine ch_get_blk_ptr &
          (blkID,gridDataStruct,dataPtr) &
          bind(c,name='ch_get_blk_ptr')
       use iso_c_binding, only : c_int, c_ptr
       integer(c_int), value, intent(in) :: blkID, gridDataStruct
       type(c_ptr), intent(out) :: dataPtr
     end subroutine ch_get_blk_ptr
  end interface

  interface
     subroutine ch_get_block_ids &
          (blkType,level,blkIDs,pNumBlks) &
          bind(c,name="ch_get_block_ids")
       use iso_c_binding, only : c_int
       integer(c_int), value, intent(in) :: blkType, level
       integer(c_int), dimension(MAXBLOCKS), intent(out) :: blkIDs
       integer(c_int), intent(out) :: pNumBlks
     end subroutine ch_get_block_ids
  end interface

  interface
     subroutine ch_get_num_blocks &
          (pNumBlks) &
          bind(c,name="ch_get_num_blocks")
       use iso_c_binding, only : c_int
       integer(c_int), intent(out) :: pNumBlks
     end subroutine ch_get_num_blocks
  end interface

  interface
     subroutine ch_get_cell_coords &
          (blkID,axis,edge,size,guardcell,coords) &
          bind(c,name="ch_get_cell_coords")
       use iso_c_binding, only : c_int, c_double
       integer(c_int), value, intent(IN) :: blkID
       integer(c_int), value, intent(IN) :: axis
       integer(c_int), value, intent(IN) :: edge
       integer(c_int), value, intent(IN) :: size
       integer(c_int), value, intent(IN) :: guardcell
       real(c_double), dimension(size), intent(OUT) :: coords
     end subroutine ch_get_cell_coords
  end interface

  interface
     subroutine ch_get_single_cell_coords &
          (blkID,edge,beginCount,ind,coords) &
          bind(c,name="ch_get_single_cell_coords")
       use iso_c_binding, only : c_int, c_double
       integer(c_int), value, intent(IN) :: blkID
       integer(c_int), value, intent(IN) :: edge
       integer(c_int), value, intent(IN) :: beginCount
       integer(c_int), dimension(MDIM), intent(IN) :: ind
       real(c_double), dimension(MDIM), intent(OUT) :: coords
     end subroutine ch_get_single_cell_coords
  end interface

  interface
     subroutine ch_get_box_info &
          (blkID,gridDataStruct,boxInfo) &
          bind(c,name="ch_get_box_info")
       use iso_c_binding, only : c_int
       use flash_ftypes, only : box_info_t
       integer(c_int), value, intent(in) :: blkID
       integer(c_int), value, intent(in) :: gridDataStruct
       type(box_info_t), intent(out) :: boxInfo
     end subroutine ch_get_box_info
  end interface

  interface
     subroutine ch_reflux &
          () &
          bind(c,name='ch_reflux')
     end subroutine ch_reflux
  end interface

  interface
     subroutine ch_post_time_step &
          () &
          bind(c,name='ch_post_time_step')
     end subroutine ch_post_time_step
  end interface

  interface
     subroutine ch_zero_flux_data &
          () &
          bind(c,name='ch_zero_flux_data')
     end subroutine ch_zero_flux_data
  end interface

  interface
     subroutine ch_put_flux_data &
          (pFluxes, dt, blkID, axis) &
          bind(c,name='ch_put_flux_data')
       use iso_c_binding, only : c_int, c_double, c_ptr
       type(c_ptr), value, intent(in) :: pFluxes
       real(c_double), value, intent(in) :: dt
       integer(c_int), value, intent(in) :: blkID
       integer(c_int), value, intent(in) :: axis
     end subroutine ch_put_flux_data
  end interface
end module chombo_f_c_interface
