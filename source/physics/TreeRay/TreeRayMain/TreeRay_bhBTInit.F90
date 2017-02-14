!!****if* source/physics/TreeRay/TreeRayMain/TreeRay_bhBTInit
!!
!! NAME
!!
!!  TreeRay_bhBTInit
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

subroutine TreeRay_bhBTInit()
  use TreeRay_data, ONLY : tr_BhUseTreeRay
  !use tr_osInterface, ONLY : tr_osBTInit
  use tr_odInterface, ONLY : tr_odBTInit
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  if (.not. tr_bhUseTreeRay) return

  !call tr_osBTInit()
  call tr_odBTInit()

  return
end subroutine TreeRay_bhBTInit

