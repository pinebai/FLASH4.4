!!****f* source/physics/Gravity/Gravity_bhNodeContrib
!!
!! NAME
!!
!!  Gravity_bhNodeContrib
!!
!!
!! SYNOPSIS
!!
!!   call Gravity_bhNodeContrib(
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
!!  Calculates contribution of the node to the gravitational
!!  potential and adds it to the solnData array (with index grv_bhGpotVar).
!!
!! ARGUMENTS
!!
!!  node        : array of the node of the tree, whose
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

subroutine Gravity_bhNodeContrib(node, trLevel, refLevel, dr, blockno, point, blkLimits, solnData)
  implicit none
#include "constants.h"
#include "Flash.h"
#include "FortranLangFeatures.fh"
  real, dimension(:), intent(IN) :: node
  integer, intent(IN) :: trLevel, refLevel, blockno
  integer, dimension(MDIM), intent(IN) :: point
  real, dimension(MDIM+2), intent(IN) :: dr
  integer, dimension(2,MDIM), intent(IN)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData

  return
end subroutine Gravity_bhNodeContrib
