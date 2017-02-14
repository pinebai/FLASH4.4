!!****if* source/physics/TreeRay/TreeRayMain/TreeRay_bhAccNode
!!
!! NAME
!!
!!  TreeRay_bhAccNode
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

subroutine TreeRay_bhAccNode(subnode, accnode)
  use TreeRay_data, ONLY : tr_BhUseTreeRay
  !use tr_osInterface, ONLY : tr_osAccNode
  use tr_odInterface, ONLY : tr_odAccNode
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  real, dimension(:), intent(IN)  :: subnode
  real, dimension(:), intent(OUT) :: accnode


  if (.not. tr_bhUseTreeRay) return

  !call tr_osAccNode(subnode, accnode)
  call tr_odAccNode(subnode, accnode)

end subroutine TreeRay_bhAccNode
