!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhBotNodeContrib
!!
!! NAME
!!
!!  gr_bhBotNodeContrib
!!
!!
!! SYNOPSIS
!!
!!  call gr_bhBotNodeContrib(
!!                          real(in)       :: node(:),
!!                          integer(in)    :: trLevel,
!!                          integer(in)    :: refLevel,
!!                          real(in)       :: dr(MDIM+2),
!!                          integer(in)    :: blockno,
!!                          integer(in)    :: point(MDIM),
!!                          integer(in)    :: blkLimits(2,MDIM),
!!                          real,pointer   :: solnData(:,:,:,:)
!!        )
!!
!! DESCRIPTION
!!
!!  Calculates contribution of the bottom-most node. Calls *_bhNodeContrib 
!!  subroutines of physical modules (recently Gravity and TreeRay) for specific 
!!  calculation. In future, it will be merged with gr_bhNodeContrib - in
!!  physical modules it has already happen in this version.
!!
!! ARGUMENTS
!!
!!  node        : array of the bottom-most (leaf) node of the tree, whose
!!                contribution is added to the gravitational potential
!!  trLevel     : level of the block-tree
!!  refLevel    : refinement level of the block with the node
!!  dr          : (1:MDIM) - position vector from the point-of-calculation to the node
!!                (MDIM+1) - square of the magnitude of the position vector
!!                (MDIM+2) - inverted magnitude of the position vector
!!  blockno     : number of block into which the point-of-calculation belongs
!!  point       : indeces of the point-of-calculation in the block
!!  blkLimits   : limits of indeces in the block
!!  solnData    : solution data from the grid
!!
!!
!!***

subroutine gr_bhBotNodeContrib(node, trLevel, refLevel, dr, blockno, point, blkLimits, solnData)
  use gr_bhData, ONLY : gr_bhOAAvg, gr_bhOAMin, gr_bhOAMax, gr_bhOACnt, &
    gr_bhTreeNodeSize
  use Gravity_interface, ONLY : Gravity_bhNodeContrib
  use TreeRay_interface, ONLY : TreeRay_bhNodeContrib
  implicit none
#include "FortranLangFeatures.fh"
#include "constants.h"
  real, dimension(:), intent(IN) :: node
  integer, intent(IN) :: trLevel, refLevel
  real, dimension(MDIM+2), intent(IN) :: dr
  integer, intent(IN) :: blockno
  integer, dimension(MDIM), intent(IN) :: point
  integer, dimension(2,MDIM), intent(IN) :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
  real :: oa

  ! call _bhBotNodeContrib subroutines of physical modules
  call Gravity_bhNodeContrib(node, trLevel, refLevel, dr, blockno, point, blkLimits, solnData)
  call TreeRay_bhNodeContrib(node, trLevel, refLevel, dr, blockno, point, blkLimits, solnData)

  ! calculate opening angle - for monitoring
  oa = gr_bhTreeNodeSize(trLevel+refLevel)*dr(MDIM+2)
  gr_bhOAAvg = gr_bhOAAvg + oa
  if (oa < gr_bhOAMin) gr_bhOAMin = oa
  if (oa > gr_bhOAMax) gr_bhOAMax = oa
  gr_bhOACnt = gr_bhOACnt + 1
  
  return
end subroutine gr_bhBotNodeContrib


