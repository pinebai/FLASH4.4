!!****if* source/physics/TreeRay/TreeRayMain/tr_bhRadToGas
!!
!! NAME
!!
!!  tr_bhRadToGas
!!
!!
!! SYNOPSIS
!!
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!
!! RESULT
!!
!!
!!***

subroutine tr_bhRadToGas(solnPoint, vol_poc, area_poc, dt)
  !use tr_osInterface, ONLY : tr_osRadToGas
  !use tr_rpInterface, ONLY : tr_rpradToGas
  !use tr_odInterface, ONLY : tr_rpradToGas
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Eos.h"
  real, DIMENSION(:), POINTER, intent(INOUT) :: solnPoint
  real, intent(IN) :: vol_poc, area_poc, dt

  !call tr_osRadToGas(solnPoint, vol_poc, area_poc)
  !call tr_rpRadToGas(solnPoint, dt)
  !call tr_odRadToGas(solnPoint, dt)

  return
end subroutine tr_bhRadToGas

