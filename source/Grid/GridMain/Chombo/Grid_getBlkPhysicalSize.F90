!!****if* source/Grid/GridMain/Chombo/Grid_getBlkPhysicalSize
!!
!! NAME
!!  Grid_getBlkPhysicalSize
!!
!! SYNOPSIS
!!
!!  Grid_getBlkPhysicalSize(integer(IN)  :: blockId,
!!                          real(OUT) :: blockSize(MDIM))
!!  
!! DESCRIPTION 
!!  Gets the physical size of the block along each dimension by calculating
!!  the deltas of each direction times the number of zones.  It is important
!!  to calculate the size of a block consistently because rounding errors
!!  can occur in some cases.  Use this function to get the size of the 
!!  block rather than calculating it in a new way.  
!!
!! ARGUMENTS
!!  blockId -local block number
!!  blockSize - returned array of size MDIM holding the size of 
!!              each dimension of the block
!!
!! EXAMPLE
!!
!!    If there is 2-d block with 8 cells along IAXIS and 4 cells along
!!    JAXIS, and a single cell size (or delta) along each dimension is 0.1, 
!!    then the part of the physical domain occupied by the block is 
!!    0.8 x 0.4 in size, as in the block shown below
!!    
!!     - - - - - - - - (1.0,0.8)
!!    | | | | | | | | |
!!     - - - - - - - - 
!!    | | | | | | | | |
!!     - - - - - - - - 
!!    | | | | | | | | |
!!     - - - - - - - - 
!!    | | | | | | | | |
!!     - - - - - - - - 
!!   (0.2,0.4)
!!  
!!
!! NOTES
!!
!!  The physical block sizes returned by this routine are
!!  Calculated using deltas.
!!
!!  This interface is general enough to work even
!!  when the block sizes (in terms of number of cells) are not uniform.
!!  
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine Grid_getBlkPhysicalSize(blockID, blockSize)
  use iso_c_binding, ONLY : c_int
  use flash_ftypes, ONLY : box_info_t
  use chombo_f_c_interface, ONLY : ch_get_box_info
  implicit none
  integer,intent(in) :: blockID
  real,dimension(MDIM),intent(out) :: blockSize
  type(box_info_t) :: boxInfo
  integer(c_int) :: blkID, gds

  blkID = blockID
  gds = CENTER
  call ch_get_box_info(blkID, gds, boxInfo);
  blockSize = boxInfo % bsize
end subroutine Grid_getBlkPhysicalSize
