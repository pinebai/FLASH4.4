!!****f* source/physics/Gravity/Gravity_bhPartErr
!!
!! NAME
!!
!!  Gravity_bhPartErr
!!
!!
!! SYNOPSIS
!!
!!  call Gravity_bhPartErr(
!!             real(in)                       :: node(:),
!!             real(in)                       :: ndSize,
!!             real(in)                       :: dr(MDIM+2),
!!             real(out)                      :: perr
!!             )
!!
!! DESCRIPTION
!!
!!  Calculates error made at a contribution of the node to the gravitational
!!  potential according to Eq. 9 in Salmon & Warren (1994, J.Comp.Phys., 136, 155).
!!  Used by the SumSquare MAC which is not fully implemented in this version.
!!
!! ARGUMENTS
!!
!!  node        : array of the node of the tree, whose
!!                contribution is added to the gravitational potential
!!  ndSize      : physical size of the node (the largest extent of the node)
!!  dr          : (1:MDIM) - position vector from the point-of-calculation to the node
!!                (MDIM+1) - square of the magnitude of the position vector
!!                (MDIM+2) - inverted magnitude of the position vector
!!  perr        : calculated error
!!
!!
!!***

subroutine Gravity_bhPartErr(node, ndSize, dr, perr)
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  real, dimension(:), intent(IN) :: node
  real, intent(IN) :: ndSize
  real, dimension(MDIM+2), intent(IN) :: dr
  real, intent(OUT) :: perr
  perr = 0.0
  return
end subroutine Gravity_bhPartErr


