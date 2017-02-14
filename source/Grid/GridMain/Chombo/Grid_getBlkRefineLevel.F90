!!****if* source/Grid/GridMain/Chombo/Grid_getBlkRefineLevel
!!
!! NAME
!!  Grid_getBlkRefineLevel
!!
!! SYNOPSIS
!!
!!  Grid_getBlkRefineLevel(integer(IN)  :: blockID,
!!                         integer(OUT) :: refineLevel)
!!  
!! DESCRIPTION 
!!  Get the refinement level of a given block as denoted by blockID.
!!  For the UG this is always 1.
!!
!! ARGUMENTS
!!  blockID - the local block number
!!  refineLevel - returned value, refinement level of block
!!
!! 
!!
!!***

#include "constants.h"

subroutine Grid_getBlkRefineLevel(blockID,refineLevel)
  use iso_c_binding, ONLY : c_int
  use flash_ftypes, ONLY : box_info_t
  use chombo_f_c_interface, ONLY : ch_get_box_info
  implicit none
  integer,intent(in) :: blockID
  integer,intent(out) :: refineLevel
  type(box_info_t) :: boxInfo
  integer(c_int) :: blkID, gds

  blkID = blockID
  gds = CENTER
  call ch_get_box_info(blkID, gds, boxInfo)
  refineLevel = boxInfo % lrefine
end subroutine Grid_getBlkRefineLevel
