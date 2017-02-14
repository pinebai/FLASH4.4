!!****f* source/physics/Gravity/Gravity_bhNormalizeNode
!!
!! NAME
!!
!!  Gravity_bhNormalizeNode
!!
!!
!! SYNOPSIS
!!
!!   call Gravity_bhNormalizeNode(
!!                          real(in)       :: smr(MDIM),
!!                          real(inout)    :: node(:)
!!        )
!!
!! DESCRIPTION
!!
!!  Called during tree build after all subnodes are added to this node.
!!  Finish calculation of the B2 coefficient, if present in the node.
!!
!! ARGUMENTS
!!
!!  smr   : sum of mass times position of subnodes
!!  node  : array of the node
!!
!!
!!***

subroutine Gravity_bhNormalizeNode(smr, node)
  implicit none
#include "constants.h"
#include "Flash.h"
  real, dimension(MDIM), intent(IN) :: smr
  real, dimension(:), intent(INOUT) :: node

  return
end subroutine Gravity_bhNormalizeNode

