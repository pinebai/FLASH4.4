!!****if* source/Grid/GridMain/UG/Grid_getBlkCornerID
!!
!! SYNOPSIS
!!
!!  call Grid_getBlkCornerID(integer(IN)  :: blockId,
!!                           integer(OUT) :: cornerID(MDIM),
!!                           integer(OUT) :: stride(MDIM),
!!                 optional, integer(OUT) :: cornerIDHigh(MDIM),
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
!!     
!!  stride  :: spacing factor between indices, in UG stride is always = 1.
!!             For PARAMESH, stride may be more than 1, depending
!!             on how far down you are in the tree.
!!
!!  cornerIDHigh :: global integer indices of the last interior zone
!!              of the block
!!
!!  inRegion :: optional argument, not valid for this implementation
!!
!! EXAMPLE
!!
!!  In a 1 dimensional UG case with 2 blocks and nxb=8
!!  The cornerID for block 1 = 1 and the cornerID for block 2 = 9 
!!  
!! 
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine Grid_getBlkCornerID(blockId, cornerID, stride, cornerIDHigh, inRegion)
  use Grid_data, ONLY : gr_blkCornerID,gr_lIndexSize, gr_axisMe
  implicit none

#include "constants.h"
#include "Flash.h"
  integer,intent(IN) :: blockId
  integer,dimension(MDIM), intent(OUT) :: cornerID, stride
  integer, optional, dimension(MDIM),intent(OUT) :: cornerIDHigh
  logical, optional, intent(IN) :: inRegion

  cornerID = gr_blkCornerID
  stride = 1
  if(present(cornerIDHigh))cornerIDHigh=(gr_axisMe+1)*gr_lIndexSize

end subroutine Grid_getBlkCornerID

