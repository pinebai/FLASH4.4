!!****f* source/physics/TreeRay/TreeRay_bhPostprocNode
!!
!! NAME
!!
!!  TreeRay_bhPostprocNode
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

subroutine TreeRay_bhPostprocNode(ndSize, node)
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  real, intent(IN) :: ndSize
  real, dimension(:), intent(INOUT)  :: node

  return
end subroutine TreeRay_bhPostprocNode

