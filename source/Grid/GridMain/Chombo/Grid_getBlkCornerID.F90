!!****if* source/Grid/GridMain/Chombo/Grid_getBlkCornerID
!!
!! SYNOPSIS
!!
!!  call Grid_getBlkCornerID(integer(IN)  :: blockId,
!!                           integer(OUT) :: cornerID(MDIM),
!!                           integer(OUT) :: stride(MDIM)
!!                  optional,integer(OUT) :: cornerIDHigh(MDIM),
!!                   optional,logical(IN) :: inRegion)
!!  
!! DESCRIPTION 
!! 
!!  Returns the global integer indices of the left most interior zone
!!  of the block and the stride of indices along each dimension.
!!  Together the cornerID and the stride make a unique identifier for
!!  each block on the grid.
!!
!!  Since Uniform Grid does not have different levels of refinement,
!!  the value of stride is always one and the cornerID can be
!!  calculated as me(axis)*blkLimits(axis), where me(axis) is the
!!  processor number along desired axis and blkLimits(axis) is the
!!  size of the blocks along the same axis.
!!
!!  CornerID counting starts at 1 and only the interior cells (no
!!  guardcells) are used to calculate the cornerID.
!! 
!! 
!! ARGUMENTS 
!!
!!  blockId :: the local blockID
!!  cornerID :: global integer indices of start of the interior zone
!!              of the block
!!  stride  :: spacing factor between indices, in UG stride is always = 1.
!!             For PARAMESH, stride may be more than 1, depending
!!             on how far down you are in the tree.
!!  cornerIDHigh :: global integer indices of the last interior zone
!!              of the block.
!!              This optional argument is IGNORED IN THIS IMPLEMENTATION.
!!  inRegion :: if present and true, the corner ID is computed relative to the region
!!              defined by the Grid scope variable gr_region
!!
!! EXAMPLE
!!
!!  In a 1 dimensional UG case with 2 blocks and nxb=8
!!  The cornerID for block 1 = 1 and the cornerID for block 2 = 9 
!!  
!! 
!!***

#include "constants.h"
#include "Flash.h"

subroutine Grid_getBlkCornerID(blockId, cornerID, stride, cornerIDHigh, inRegion)
  use iso_c_binding, ONLY : c_int
  use flash_ftypes, ONLY : box_info_t
  use chombo_f_c_interface, ONLY : ch_get_box_info
  implicit none
  integer,intent(IN) :: blockId
  integer,dimension(MDIM), intent(OUT) :: cornerID, stride
  integer,dimension(MDIM),optional,intent(OUT) :: cornerIDHigh
  logical, optional, intent(IN) :: inRegion
  type(box_info_t) :: boxInfo
  integer(c_int) :: blkID, gds

  blkID = blockId
  gds = CENTER
  call ch_get_box_info(blkID, gds, boxInfo)
  cornerID = boxInfo % corner
  stride = boxInfo % stride
  !! DEV there is no implementation to provide value into optional arguement cornerIDHigh
end subroutine Grid_getBlkCornerID
