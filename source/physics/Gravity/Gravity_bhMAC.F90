!!****f* source/physics/Gravity/Gravity_bhMAC
!!
!! NAME
!!
!!  Gravity_bhMAC
!!
!!
!! SYNOPSIS
!!
!!   logical res = Gravity_bhMAC(
!!                          real(in)       :: node(:),
!!                          real(in)       :: ndSize2,
!!                          real(in)       :: dr(MDIM+2),
!!                          integer(in)    :: blockno,
!!                          integer(in)    :: point(MDIM),
!!                          integer(in)    :: blkLimits(2,MDIM),
!!                          real,pointer   :: solnData(:,:,:,:)
!!        )
!!
!! DESCRIPTION
!!
!!  Multipole Acceptance Criterion. Determines whether the contribution of the
!!  node to the potential at the point of calculation will have a nacessary
!!  accuracy.
!!
!! ARGUMENTS
!!
!!  node        : array of the node tested
!!  ndSize2     : square of the physical size of the node (the largest 
!!                extent of the node)
!!  dr          : (1:MDIM) - position vector from the point-of-calculation to the node
!!                (MDIM+1) - square of the magnitude of the position vector
!!                (MDIM+2) - inverted magnitude of the position vector
!!  blockno     : number of block into which the point-of-calculation belongs
!!  point       : indeces of the point-of-calculation in the block
!!  blkLimits   : limits of indeces in the block
!!  solnData    : solution data from the grid
!!
!! RESULT
!!  Returns TRUE if the node is accepted for the calculation of the potential.
!!  Otherwise returns FALSE.
!!
!!
!!***

logical function Gravity_bhMAC(node, ndSize2, dr, blockno, point, blkLimits, solnData)
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  real, dimension(:), intent(IN) :: node
  real, intent(IN) :: ndSize2
  integer, dimension(MDIM), intent(IN) :: point
  real, dimension(MDIM+2), intent(IN) :: dr
  integer, intent(IN) :: blockno
  integer, dimension(2,MDIM), intent(IN)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER :: solnData

  Gravity_bhMAC = .true.
  return
end function Gravity_bhMAC


