!!****if* source/physics/TreeRay/TreeRayMain/tr_odFinalizeCell
!!
!! NAME
!!
!!  tr_odFinalizeCell
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

subroutine tr_odFinalizeCell(solnpoint, vol_poc, cdMaps)
  use TreeRay_data, ONLY : tr_nCd, tr_nPix
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  real, DIMENSION(:), POINTER, intent(INOUT) :: solnPoint
  real, intent(IN) :: vol_poc
  real, DIMENSION(tr_nCd, 0:tr_nPix-1), intent(IN) :: cdMaps

  return
end subroutine tr_odFinalizeCell

