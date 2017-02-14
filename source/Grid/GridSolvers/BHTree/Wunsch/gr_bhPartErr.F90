!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhPartErr
!!
!! NAME
!!
!!  gr_bhPartErr
!!
!!
!! SYNOPSIS
!!
!!  call gr_bhPartErr(
!!             real(in)    :: node(:),
!!             real(in)    :: ndSize,
!!             real(in)    :: dr(MDIM+2),
!!             real(out)   :: perr(GR_TREE_NMODULES)
!!             )
!!
!! DESCRIPTION
!!
!!  Calculates errors made at contributions of the node for individual physical
!!  modules by calling corresponding *_bhPartErr subroutines.
!!  Resulting errors are collected in the perr array.
!!
!! ARGUMENTS
!!
!!  node        : array of the node of the tree, whose
!!                contribution is added to the gravitational potential
!!  ndSize      : physical size of the node (the largest extent of the node)
!!  dr          : (1:MDIM) - position vector from the point-of-calculation to the node
!!                (MDIM+1) - square of the magnitude of the position vector
!!                (MDIM+2) - inverted magnitude of the position vector
!!  perr        : calculated errors
!!
!!
!!***

subroutine gr_bhPartErr(node, ndSize, dr, perr)
  use gr_bhData, ONLY : gr_bhUseGravity, gr_bhUseTreeRay, &
    GR_TREE_NMODULES, GR_TREE_IGRAVITY, GR_TREE_ITREERAY
  use Gravity_interface, ONLY : Gravity_bhPartErr
  use TreeRay_interface, ONLY : TreeRay_bhPartErr
  implicit none
#include "constants.h"
  real, dimension(:), intent(IN) :: node
  real, intent(IN) :: ndSize
  real, dimension(MDIM+2), intent(IN) :: dr
  real, dimension(GR_TREE_NMODULES), intent(OUT) :: perr

  perr = 0.0

  ! call _bhPartErr subroutines of physical modules
  if (gr_bhUseGravity) &
  & call Gravity_bhPartErr(node, ndSize, dr, perr(GR_TREE_IGRAVITY))
  if (gr_bhUseTreeRay) &
  & call TreeRay_bhPartErr(node, ndSize, dr, perr(GR_TREE_ITREERAY))

  return
end subroutine gr_bhPartErr

