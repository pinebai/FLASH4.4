!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhAccNode
!!
!! NAME
!!
!!  gr_bhAccNode
!!
!!
!! SYNOPSIS
!!
!!  call gr_bhAccNode(
!!                          real(in)       :: subnode(:),
!!                          real(out)      :: accnode(:)
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

subroutine gr_bhAccNode(subnode, accnode)
  use gr_bhData, ONLY : GR_TREE_IM, GR_TREE_IX, GR_TREE_IY, GR_TREE_IZ
  use Gravity_interface, ONLY : Gravity_bhAccNode
  use TreeRay_interface, ONLY : TreeRay_bhAccNode
  implicit none
  real, dimension(:), intent(IN)  :: subnode
  real, dimension(:), intent(INOUT) :: accnode

  ! adds mass and mass*position contributions of subnode to accnode
  accnode(GR_TREE_IM) = accnode(GR_TREE_IM) + subnode(GR_TREE_IM)
  accnode(GR_TREE_IX) = accnode(GR_TREE_IX) + subnode(GR_TREE_IM)*subnode(GR_TREE_IX)
  accnode(GR_TREE_IY) = accnode(GR_TREE_IY) + subnode(GR_TREE_IM)*subnode(GR_TREE_IY)
  accnode(GR_TREE_IZ) = accnode(GR_TREE_IZ) + subnode(GR_TREE_IM)*subnode(GR_TREE_IZ)

  ! call _bhAccNode subroutines of physical modules
  call Gravity_bhAccNode(subnode, accnode)
  call TreeRay_bhAccNode(subnode, accnode)

  return
end subroutine gr_bhAccNode
