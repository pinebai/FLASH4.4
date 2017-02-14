!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhSelfContrib
!!
!! NAME
!!
!!  gr_bhSelfContrib
!!
!!
!! SYNOPSIS
!!
!!  call gr_bhSelfContrib(
!!                  real(in)       :: node(:),
!!                  real(in)       :: cellsize(MDIM),
!!                  integer(in)    :: blockno,
!!                  integer(in)    :: point(MDIM),
!!                  integer(in)    :: blkLimits(2,MDIM),
!!                  real,pointer   :: solnData(:,:,:,:)
!!        )
!!
!! DESCRIPTION
!!
!!  Calculates contribution of the cell at the point-of-calculation.
!!
!! ARGUMENTS
!!
!!  node        : array of the node of the tree, whose
!!                contribution is added to the gravitational potential
!!  cellsize    : physical size of the cell (in each dimension)
!!  blockno     : number of block into which the point-of-calculation belongs
!!  point       : indeces of the point-of-calculation in the block
!!  blkLimits   : limits of indeces in the block
!!  solnData    : solution data from the grid
!!
!!***

subroutine gr_bhSelfContrib(node, cellsize, blockno, point, blkLimits, solnData)
  use Gravity_interface, ONLY : Gravity_bhSelfContrib
  use TreeRay_interface, ONLY : TreeRay_bhSelfContrib
  implicit none
#include "FortranLangFeatures.fh"
#include "constants.h"
  real, dimension(:), intent(IN) :: node
  real, dimension(MDIM), intent(IN) :: cellsize
  integer, intent(IN) :: blockno
  integer, dimension(MDIM), intent(IN) :: point
  integer, dimension(2,MDIM), intent(IN) :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData

  ! call _bhSelfContrib subroutines of physical modules
  call Gravity_bhSelfContrib(node, cellsize, blockno, point, blkLimits, solnData)
  call TreeRay_bhSelfContrib(node, cellsize, blockno, point, blkLimits, solnData)
  
  return
end subroutine gr_bhSelfContrib


