!!****f* source/physics/TreeRay/TreeRay_bhPartErr
!!
!! NAME
!!
!!  TreeRay_bhPartErr
!!
!!
!! SYNOPSIS
!!
!!  call TreeRay_bhPartErr(
!!             real, dimension(:), intent(IN)                 :: node
!!             real, intent(IN)                               :: ndSize  
!!             real, dimension(MDIM+2), intent(IN)            :: dr
!!             real, intent(OUT)                              :: perr
!!             )
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!  node
!!  ndSize  
!!  dr
!!  blockno
!!  point
!!
!! RESULT
!!
!!
!!***

subroutine TreeRay_bhPartErr(node, ndSize, dr, perr)
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  real, dimension(:), intent(IN) :: node
  real, intent(IN) :: ndSize
  real, dimension(MDIM+2), intent(IN) :: dr
  real, intent(OUT) :: perr
  perr = 0.0
  return
end subroutine TreeRay_bhPartErr


