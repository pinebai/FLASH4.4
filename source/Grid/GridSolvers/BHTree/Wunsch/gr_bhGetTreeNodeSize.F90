!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhGetTreeNodeSize
!!
!! NAME
!!
!!  gr_bhGetTreeNodeSize
!!
!! 
!! SYNOPSIS
!!
!!  real nSize = call gr_bhGetTreeNodeSize(
!!                          integer(in)    :: level
!!            )
!!
!!
!! DESCRIPTION
!!
!!  Calculates size (along the longest edge) of a tree node at a given level.
!!
!! ARGUMENTS
!!
!!  level     : level of the node, 1 = refinement level 1
!!
!! RESULT
!!
!!  Returns physical size of the tree node - length of the longest edge.
!!
!!***

real function gr_bhGetTreeNodeSize(level)
  use gr_bhData, ONLY : gr_bhTreeLrefineMax, gr_bhTreeLevels
  use Grid_interface, ONLY : Grid_getMinCellSizes

  implicit none
#include "Flash.h"
#include "constants.h"
  integer, intent(IN) :: level
  real :: mcs(MDIM), max_mcs

  call Grid_getMinCellSizes(mcs)
  max_mcs = max(mcs(IAXIS), MCS(JAXIS), mcs(KAXIS))
  gr_bhGetTreeNodeSize = max_mcs * 2**(gr_bhTreeLrefineMax+gr_bhTreeLevels-level)

  return
end function gr_bhGetTreeNodeSize
