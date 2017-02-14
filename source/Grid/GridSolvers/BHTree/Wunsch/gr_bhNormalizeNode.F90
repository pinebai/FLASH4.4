!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhNormalizeNode
!!
!! NAME
!!
!!  gr_bhNormalizeNode
!!
!!
!! SYNOPSIS
!!
!!  call gr_bhNormalizeNode(
!!                          real(inout)    :: node(:)
!!        )
!!
!! DESCRIPTION
!!
!!  Called during tree build after all subnodes are added to this node.
!!  Calculates mass centre of the node and calls *_bhNormalizeNode subroutines
!!  of physical modules.
!!
!! ARGUMENTS
!!
!!  node  : array of the node
!!
!!
!!***

subroutine gr_bhNormalizeNode(node)
  use gr_bhData, ONLY : GR_TREE_IM, GR_TREE_IX, GR_TREE_IY, GR_TREE_IZ
  use Gravity_interface, ONLY : Gravity_bhNormalizeNode
  use TreeRay_interface, ONLY : TreeRay_bhNormalizeNode
  implicit none
#include "constants.h"
  real, dimension(:), intent(INOUT) :: node
  real, dimension(MDIM) :: smr

  ! store sum of m_i * r{i,mc} for subroutines of physical modules
  smr(IAXIS) = node(GR_TREE_IX)
  smr(JAXIS) = node(GR_TREE_IY)
  smr(KAXIS) = node(GR_TREE_IZ)

  ! calculate mass centre position of the node
  node(GR_TREE_IX) = node(GR_TREE_IX) / node(GR_TREE_IM)
  node(GR_TREE_IY) = node(GR_TREE_IY) / node(GR_TREE_IM)
  node(GR_TREE_IZ) = node(GR_TREE_IZ) / node(GR_TREE_IM)

  ! call _bhNormalizeNode subroutines of physical modules
  call Gravity_bhNormalizeNode(smr, node)
  call TreeRay_bhNormalizeNode(smr, node)

  return
end subroutine gr_bhNormalizeNode

