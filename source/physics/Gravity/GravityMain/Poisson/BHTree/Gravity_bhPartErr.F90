!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/Gravity_bhPartErr
!!
!! NAME
!!
!!  Gravity_bhPartErr
!!
!!
!! SYNOPSIS
!!
!!  call Gravity_bhPartErr(
!!             real(in)    :: node(:),
!!             real(in)    :: ndSize,
!!             real(in)    :: dr(MDIM+2),
!!             real(out)   :: perr
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
  use Gravity_data, ONLY : useGravity, grv_bhIB2, grv_bhIB3, grv_bhAccErr, &
    grv_bhNewton, grv_bhIM
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  real, dimension(:), intent(IN) :: node
  real, intent(IN) :: ndSize
  real, dimension(MDIM+2), intent(IN) :: dr
  real, intent(OUT) :: perr
  real :: ar ! node aspect ratio
  real :: B3_low ! low bound on the B3 moment

  if (.not. useGravity) return

  !print *, "grv_bhPE: ", node, "|", ndSize, "|", dr
  ar = ndSize*dr(MDIM+2)
  if (ar < 1.0) then
    B3_low = sqrt(node(grv_bhIB2)**3/node(grv_bhIM))
    perr = 1.0 / (dr(MDIM+1)*(1.0-ar)*(1.0-ar)) &
    &    * grv_bhNewton * (3.0*node(grv_bhIB2)*dr(MDIM+2)*dr(MDIM+2) &
    &    -  2.0*B3_low*dr(MDIM+2)*dr(MDIM+2)*dr(MDIM+2))
    if (perr /= perr) then
      print *, "PERR NAN: ", node, "|", dr, "|", ar, ndSize
    endif
  else
    perr = 2*grv_bhAccErr
  endif
  perr = min(2*grv_bhAccErr, abs(perr))
  !print *, "PERR: ", perr, ar

  return
end subroutine Gravity_bhPartErr


