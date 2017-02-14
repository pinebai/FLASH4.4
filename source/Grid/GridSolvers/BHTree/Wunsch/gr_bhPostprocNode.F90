!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhPostprocNode
!!
!! NAME
!!
!!  gr_bhPostprocNode
!!
!!
!! SYNOPSIS
!!
!!  call gr_bhPostprocNode(
!!             real(in)                       :: ndSize,
!!             real(inout)                    :: node(:)
!!             )
!!
!! DESCRIPTION
!!
!!  Called after the tree is built for each node.
!!
!! ARGUMENTS
!!
!!  ndSize      : physical size of the node (the largest extent of the node)
!!  node        : array of the tree node which is processed
!!
!!
!!***

subroutine gr_bhPostprocNode(ndSize, node)
  use Gravity_interface, ONLY : Gravity_bhPostprocNode
  use TreeRay_interface, ONLY : TreeRay_bhPostprocNode
  implicit none
  real, intent(IN) :: ndSize
  real, dimension(:), intent(INOUT) :: node

  ! call _bhPostprocNode subroutines of physical modules
  call Gravity_bhPostprocNode(ndSize, node)
  call TreeRay_bhPostprocNode(ndSize, node)

  return
end subroutine gr_bhPostprocNode

