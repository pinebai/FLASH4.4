!!****f* source/physics/Gravity/Gravity_bhAccNode
!!
!! NAME
!!
!!  Gravity_bhAccNode
!!
!!
!! SYNOPSIS
!!
!!   call Gravity_bhAccNode(
!!                          real(in)       :: subnode(:),
!!                          real(inout)    :: accnode(:)
!!        )
!!
!! DESCRIPTION
!!
!!  Called during tree build. Adds values of subnode into accnode.
!!
!! ARGUMENTS
!!
!!  subnode     : array of the node of the tree which is added to accnode
!!  accnode     : array of the node into which subnode contribution is added
!!
!!
!!
!!***

subroutine Gravity_bhAccNode(subnode, accnode)
  implicit none
#include "constants.h"
#include "Flash.h"
  real, dimension(:), intent(IN)  :: subnode
  real, dimension(:), intent(INOUT) :: accnode
  return
end subroutine Gravity_bhAccNode
