!!****f* source/physics/TreeRay/TreeRay_bhAccNode
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
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  real, dimension(:), intent(IN)  :: subnode
  real, dimension(:), intent(INOUT) :: accnode
  return
end subroutine TreeRay_bhAccNode
