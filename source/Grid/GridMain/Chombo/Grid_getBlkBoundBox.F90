!!****if* source/Grid/GridMain/Chombo/Grid_getBlkBoundBox
!!
!! NAME
!!  Grid_getBlkBoundBox
!!
!! SYNOPSIS
!!
!! 
!!  Grid_getBlkBoundBox(integer(IN)  :: blockId,
!!                         real(OUT) :: boundBox(2, MDIM))
!!  
!! DESCRIPTION 
!!
!!  Gets the physical domain bounding box of the block identified 
!!  by blockId.  For each dimension the left (lower or forward) 
!!  physical coordinate of the block edge and the right (upper or back) 
!!  physical coordinate of the block edge is returned.  See arguments
!!  below for more detail.
!!
!! ARGUMENTS
!!
!!  blockId -local block number
!!
!!  boundBox - returned array holding the boundBox coordinates in
!!             each dimension
!!
!!            for readability, in constants.h we define IAXIS = 1, JAXIS = 2, KAXIS = 3
!!
!!            boundBox(1,IAXIS) = left edge coordinate of block in x direction
!!            boundBox(2,IAXIS) = right edge coordinate of block in x direction
!!            boundBox(1,JAXIS) = top edge coordinate of block in y direction
!!            boundBox(2,JAXIS) = bottom edge coordinate of block in y direction
!!            boundBox(1,KAXIS) = front edge coordinate of block in z direction
!!            boundBox(2,KAXIS) = back edge coordinate of block in z direction
!!
!! EXAMPLE
!!  
!!   In 2 dimensions, if physical coordinates are ...
!!    
!!     ________________(0.5 1.0)
!!    |                |
!!    |                |
!!    |                |
!!    |                |
!!    |                |
!!    |                |
!!    |                |
!!    |_______________ |
!!  (-0.5, 0.0)
!!
!!
!!
!!     boundBox(1, IAXIS) = -0.5
!!     boundBox(2, IAXIS) = 0.5
!!     boundBox(1, JAXIS) = 0.0
!!     boundBox(2, JAXIS) = 1.0
!!     boundBox(1, KAXIS) = 1 !returned as 1 because only 2 dims
!!     boundBox(1, KAXIS) = 1 !returned as 1 because only 2 dims
!!
!!
!!***
#include "constants.h"
#include "Flash.h"

subroutine Grid_getBlkBoundBox(blockId,boundBox)
  use iso_c_binding, ONLY : c_int
  use flash_ftypes, ONLY : box_info_t
  use chombo_f_c_interface, ONLY : ch_get_box_info
  implicit none
  integer,intent(in) :: blockId
  real,dimension(2,MDIM),intent(out) :: boundBox
  type(box_info_t) :: boxInfo
  integer(c_int) :: blkID, gds

  blkID = blockId
  gds = CENTER
  call ch_get_box_info(blkID, gds, boxInfo);
  boundBox(1,:) = boxInfo % lowBndBox(:)
  boundBox(2,:) = boxInfo % highBndBox(:)
end subroutine Grid_getBlkBoundBox
