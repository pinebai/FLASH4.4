!!****f* source/physics/Gravity/Gravity_bhFillBotNode
!!
!! NAME
!!
!!  Gravity_bhFillBotNode
!!
!!
!! SYNOPSIS
!!
!!   call Gravity_bhFillBotNode(
!!                          integer(in)    :: blockno,
!!                          integer(in)    :: point(MDIM),
!!                          integer(in)    :: blkLimits(2,MDIM),
!!                          real,pointer   :: solnData(:,:,:,:),
!!                          real(inout)    :: botnode(:)
!!        )
!!
!! DESCRIPTION
!!
!!  Called during tree build. Fills the bottom-most (leaf) node 
!!  (corresponding to a single grid cell) with values from the grid.
!!  Recently empty, because evrything is done in gr_bhFillBotNode.
!!
!! ARGUMENTS
!!
!!  blockno     : number of block into which the node belongs
!!  point       : indeces of the cell (botnode) in the block
!!  blkLimits   : limits of indeces in the block
!!  solnData    : solution data from the grid
!!  botnode     : array of the bottom-most (leaf) node of the tree
!!                which is filled by values from the grid
!!
!!***

subroutine Gravity_bhFillBotNode(blockno, point, blkLimits, solnData, botnode)
  implicit none
#include "constants.h"
  integer, intent(IN) :: blockno
  integer, dimension(MDIM), intent(IN) :: point
  integer, dimension(2,MDIM), intent(IN)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER :: solnData
  real, dimension(:), intent(INOUT) :: botnode

  return
end subroutine Gravity_bhFillBotNode
