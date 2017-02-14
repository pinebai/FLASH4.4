!!****if* source/physics/TreeRay/TreeRayMain/tr_bhFinalizeCell
!!
!! NAME
!!
!!  tr_bhFinalizeCell
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

subroutine tr_bhFinalizeCell(solnPoint, dl_poc, eflux, cdMaps)
  use TreeRay_data, ONLY : tr_nEb, tr_nPix, tr_nCd
  use tr_odInterface, ONLY : tr_odFinalizeCell
  !use tr_osInterface, ONLY : tr_osFinalizeCell
  !use tr_rpInterface, ONLY : tr_rpFinalizeCell
  implicit none
#include "constants.h"
#include "Flash.h"
  real, DIMENSION(:), POINTER, intent(INOUT) :: solnPoint
  real, DIMENSION(MDIM), intent(IN) :: dl_poc
  real, DIMENSION(tr_nEb, 0:tr_nPix-1), intent(IN) :: eflux
  real, DIMENSION(tr_nCd, 0:tr_nPix-1), intent(IN) :: cdMaps

  real, DIMENSION(tr_nEb) :: phFluxInt
  integer :: ipix, ieb
  real :: vol_poc


  vol_poc = dl_poc(IAXIS)*dl_poc(JAXIS)*dl_poc(KAXIS)
  call tr_odFinalizeCell(solnPoint, dl_poc, cdMaps)
  !call tr_osFinalizeCell(solnPoint, vol_poc, eflux, phFluxInt)
  !call tr_rpFinalizeCell(solnPoint, vol_poc, eflux, phFluxInt)

  return
end subroutine tr_bhFinalizeCell

