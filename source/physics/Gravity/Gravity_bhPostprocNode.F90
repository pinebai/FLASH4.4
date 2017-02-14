!!****f* source/physics/Gravity/Gravity_bhPostprocNode
!!
!! NAME
!!
!!  Gravity_bhPostprocNode
!!
!!
!! SYNOPSIS
!!
!!  call Gravity_bhPostprocNode(
!!             real(in)                       :: ndSize,
!!             real(inout)                    :: node(:)
!!             )
!!
!! DESCRIPTION
!!
!!  Called after the tree is built for each node. With certain MACs,
!!  calculates the minimum distance needed for this node to be accepted.
!!
!! ARGUMENTS
!!
!!  ndSize      : physical size of the node (the largest extent of the node)
!!  node        : array of the tree node which is processed
!!
!!
!!***

subroutine Gravity_bhPostprocNode(ndSize, node)
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  real, intent(IN) :: ndSize
  real, dimension(:), intent(INOUT)  :: node

  return
end subroutine Gravity_bhPostprocNode

