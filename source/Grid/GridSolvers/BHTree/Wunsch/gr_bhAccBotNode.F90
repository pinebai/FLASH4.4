!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhAccBotNode
!!
!! NAME
!!
!!  gr_bhAccBotNode
!!
!!
!! SYNOPSIS
!!
!!  call gr_bhAccBotNode(
!!                          integer(in)    :: blockno,
!!                          integer(in)    :: point(MDIM),
!!                          integer(in)    :: blkLimits(2,MDIM),
!!                          real(in)       :: locCoords(:,:,:),
!!                          real,pointer   :: solnData(:,:,:,:),
!!                          real(in)       :: botnode(:),
!!                          real(out)      :: accnode(:)
!!             )
!!
!! DESCRIPTION
!!
!!  Called during tree build. Adds values of botnode into accnode.
!!
!! ARGUMENTS
!!
!!  blockno     : number of block into which the node belongs
!!  point       : indeces of the cell (botnode) in the block
!!  blkLimits   : limits of indeces in the block
!!  locCoords   : array with coordinates of allblocks on a given cpu
!!  solnData    : solution data from the grid
!!  botnode     : array of the bottom-most (leaf) node of the tree, i.e. a grid cell
!!  accnode     : array of the node into which botnode contribution is added
!!
!!
!!***

subroutine gr_bhAccBotNode(blockno, point, blkLimits, locCoords, solnData, botnode, accnode)
  use gr_bhData, ONLY : gr_bhLocCoords, GR_TREE_IM, GR_TREE_IX, GR_TREE_IY, &
    & GR_TREE_IZ
  use Gravity_interface, ONLY : Gravity_bhAccBotNode
  use TreeRay_interface, ONLY : TreeRay_bhAccBotNode
  implicit none
#include "constants.h"
  integer, intent(IN) :: blockno
  integer, dimension(MDIM), intent(IN) :: point
  integer, dimension(2,MDIM), intent(IN)   :: blkLimits
  real, dimension(:,:,:), intent(IN) :: locCoords
  real, DIMENSION(:,:,:,:), POINTER :: solnData
  real, dimension(:), intent(IN) :: botnode
  real, dimension(:), intent(INOUT) :: accnode

  ! adds mass and mass*position contributions of botnode to accnode
  accnode(GR_TREE_IM) = accnode(GR_TREE_IM) + botnode(GR_TREE_IM)
  accnode(GR_TREE_IX) = accnode(GR_TREE_IX) + botnode(GR_TREE_IM) &
  & * gr_bhLocCoords(point(IAXIS)-blkLimits(LOW,IAXIS)+1, IAXIS, blockno)
  accnode(GR_TREE_IY) = accnode(GR_TREE_IY) + botnode(GR_TREE_IM) &
  & * gr_bhLocCoords(point(JAXIS)-blkLimits(LOW,JAXIS)+1, JAXIS, blockno)
  accnode(GR_TREE_IZ) = accnode(GR_TREE_IZ) + botnode(GR_TREE_IM) &
  & * gr_bhLocCoords(point(KAXIS)-blkLimits(LOW,KAXIS)+1, KAXIS, blockno)

  ! call _bhAccBotNode subroutines of physical modules
  call Gravity_bhAccBotNode(blockno, point, blkLimits, locCoords, solnData, botnode, accnode)
  call TreeRay_bhAccBotNode(blockno, point, blkLimits, locCoords, solnData, botnode, accnode)

  return
end subroutine gr_bhAccBotNode
