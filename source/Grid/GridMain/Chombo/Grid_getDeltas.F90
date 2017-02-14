!!****if* source/Grid/GridMain/Chombo/Grid_getDeltas
!!
!! NAME
!!  Grid_getDeltas
!!
!! SYNOPSIS
!!
!!  call Grid_getDeltas(integer(IN):: blockId, 
!!                      real(OUT)::del(MDIM))
!!  
!! DESCRIPTION 
!!  
!!  Gets the dx/dy/dz for a given blockId on the Grid
!!  dx is the size of one zone in the x direction of a block
!!
!!  
!! ARGUMENTS 
!!
!!  blockId - local block number
!!  del - array of size MDIM returned holding the dx, dy, and dz values
!!
!!  
!!***

#include "constants.h"

subroutine Grid_getDeltas(blockID,del)
  use iso_c_binding, ONLY : c_int
  use flash_ftypes, ONLY : box_info_t
  use chombo_f_c_interface, ONLY : ch_get_box_info
  implicit none
  integer, intent(IN) :: blockID
  real, dimension(MDIM),intent(OUT) :: del
  type(box_info_t) :: boxInfo
  integer(c_int) :: blkID, gds

  blkID = blockID
  gds = CENTER
  call ch_get_box_info(blkID, gds, boxInfo)
  del = boxInfo % deltas
end subroutine Grid_getDeltas
