!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhFillBotNode
!!
!! NAME
!!
!!  gr_bhFillBotNode
!!
!!
!! SYNOPSIS
!!
!!  call gr_bhFillBotNode(
!!                          integer(in)    :: blockno,
!!                          integer(in)    :: point(MDIM),
!!                          integer(in)    :: blkLimits(2,MDIM),
!!                          real,pointer   :: solnData(:,:,:,:),
!!                          real(in)       :: botnode(:)
!!        )
!!
!! DESCRIPTION
!!
!!  Called during tree build. Fills the bottom-most (leaf) node 
!!  (corresponding to a single grid cell) with values from the grid.
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

subroutine gr_bhFillBotNode(blockno, point, blkLimits, solnData, botnode)
  use Grid_interface, ONLY : Grid_getSingleCellVol
  use gr_bhData, ONLY : GR_TREE_IM, gr_bhDensVar
  use Gravity_interface, ONLY : Gravity_bhFillBotNode
  use TreeRay_interface, ONLY : TreeRay_bhFillBotNode
  implicit none
#include "FortranLangFeatures.fh"
#include "constants.h"
  integer, intent(IN) :: blockno
  integer, dimension(MDIM), intent(IN) :: point
  integer, dimension(2,MDIM)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
  real, dimension(:), intent(OUT) :: botnode
  real dvol

  ! calculates mass of the node and store it in the botnode array
  call Grid_getSingleCellVol(blockno, INTERIOR, point, dvol)
  botnode(GR_TREE_IM) = solnData(gr_bhDensVar, point(IAXIS), point(JAXIS), point(KAXIS))*dvol

  ! call _bhFillBotNode subroutines of physical modules
  call Gravity_bhFillBotNode(blockno, point, blkLimits, solnData, botnode)
  call TreeRay_bhFillBotNode(blockno, point, blkLimits, solnData, botnode)

  return
end subroutine gr_bhFillBotNode

