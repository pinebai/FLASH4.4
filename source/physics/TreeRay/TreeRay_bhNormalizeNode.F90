!!****f* source/physics/TreeRay/TreeRay_bhNormalizeNode
!!
!! NAME
!!
!!  TreeRay_bhNormalizeNode
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

subroutine TreeRay_bhNormalizeNode(smr, node)
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  real, dimension(MDIM), intent(IN) :: smr
  real, dimension(:), intent(INOUT) :: node

  return
end subroutine TreeRay_bhNormalizeNode

