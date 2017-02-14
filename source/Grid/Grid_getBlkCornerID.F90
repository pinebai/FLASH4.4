!!****f* source/Grid/Grid_getBlkCornerID
!!
!! NAME
!!  Grid_getBlkCornerID
!!
!! SYNOPSIS
!!
!!  call Grid_getBlkCornerID(integer(IN)  :: blockId,
!!                           integer(OUT) :: cornerID(MDIM),
!!                           integer(OUT) :: stride(MDIM),
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
!!  The block's cornerID is determined by calculating the global
!!  integer index of each cell on the global domain as if the grid
!!  were fully refined.  (The uniform grid is like a fully refined
!!  grid when calculating the cornerID.)  Another way to put it is
!!  that the cornerID is just the left most interior cell global index
!!  of a block in a dimension, as if the grid were fully refined.
!! 
!!  stride is defined as = 2**(lrefine_max - lrefine(blockId)) and
!!  indicates the spacing factor between one cell and one directly to
!!  its right when calculating the cornerID.
!!
!!  This means that if a block is at the highest refinement level,
!!  then it will always have a stride = 1 and the cornerID of the
!!  block to the current block's right will always be cornerID plus
!!  the index size of the block(NXB with fixed block sizes).  (The
!!  uniform grid always has a stride = 1.)
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
!!  stride  :: spacing factor between indices. In UG, stride is always = 1.
!!             For PARAMESH, stride may be more than 1, depending
!!             on how far down you are in the tree.
!!  cornerIDHigh :: global integer indices of the last interior zone
!!              of the block
!!
!!  inRegion :: if present and true, the corner ID is computed relative to the region
!!              defined by the Grid scope variable gr_region in presence of AMR
!!
!! EXAMPLE
!!
!! Example 1:
!!  In a 1 dimensional UG case with 2 blocks and nxb=8
!!  The cornerID for block 1 = 1 and the cornerID for block 2 = 9 
!!
!!
!! Example 2:
!!  In a 1 dimensional PARAMESH case with lrefine_max = 4, nxb=8, the
!!  global indices range from 1:64.  If the grid is fully refined then
!!  there are 8 blocks in the x direction and their cornerIDs are
!!  1,9,17,25,33,41,49,57 respectively. stride=1
!!  
!!  If the entire grid is at a refinement level = 3 then there are 4
!!  blocks in the x direction and their cornerIDs are 1,17,33,49
!!  respectively.  stride = 2
!!
!!  If the entire grid is at a refinement level = 2 then there are 2
!!  blocks in the x direction and their cornerIDs are 1,33
!!  respectively.  stride = 4 (meaning it takes stride*nxb to get to
!!  the next block's cornerID)
!!
!!  And so on.
!!
!!  Multiple blocks can have the same cornerID as in the above
!!  example, but they can not have the same cornerID AND stride.
!! 
!!***

subroutine Grid_getBlkCornerID(blockId, cornerID, stride,cornerIDHigh, inRegion)

  implicit none

#include "constants.h"

  integer,intent(IN)  :: blockId
  integer,dimension(MDIM), intent(OUT) :: cornerID, stride
  integer,dimension(MDIM),optional,intent(OUT) :: cornerIDHigh
  logical, optional, intent(IN) :: inRegion
  cornerID=1
  stride=1
  if(present(cornerIDHigh))cornerIDHigh=1
end subroutine Grid_getBlkCornerID
