!!****if* source/physics/TreeRay/TreeRayMain/TreeRay_bhInitFieldVar
!!
!! NAME
!!
!!  TreeRay_bhInitFieldVar
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

subroutine TreeRay_bhInitFieldVar
  use TreeRay_data, ONLY : tr_bhLocRelErr, tr_bhMaxRelEradErr, &
    tr_bhLocEradTot, tr_bhLocMionTot, tr_bhUseTreeRay
  !use tr_osInterface, ONLY : tr_osInitFieldVar
  use tr_odInterface, ONLY : tr_odInitFieldVar
  !use tr_rpInterface, ONLY : tr_rpInitFieldVar
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"


  if (.not. tr_bhUseTreeRay) return

  ! reset variables to store error - for iterations
  tr_bhLocRelErr = 0.0
  tr_bhMaxRelEradErr = 0.0
  tr_bhLocEradTot = 0.0
  tr_bhLocMionTot = 0.0


    !call tr_osInitFieldVar()
    call tr_odInitFieldVar()
    !call tr_rpInitFieldVar()


  return
end subroutine TreeRay_bhInitFieldVar
