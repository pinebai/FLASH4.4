!!****f* source/Grid/Grid_getGlobalIndexLimits
!!
!! NAME
!!  Grid_getGlobalIndexLimits
!!
!! SYNOPSIS
!!
!!  Grid_getGlobalIndexLimits(integer(OUT) :: globalIndexLimits(MDIM))
!!  
!! DESCRIPTION 
!!  Gets the integer index dimensions of the entire grid 
!!  across all processors. Guardcells are not included.
!!
!!  globalIndexLimits(IAXIS) = highest index of grid in i dir  
!!  globalIndexLimits(JAXIS) = highest index of grid in j dir   
!!  globalIndexLimits(KAXIS) = highest index of grid in k dir   
!!
!!  (IAXIS, JAXIS and KAXIS are defined in constants.h
!!  and are set to 1,2 and 3 respectively)
!!
!!  In an adaptive mesh, the highest index is returned as if
!!  the entire mesh was fully refined.
!!
!! ARGUMENTS
!!  globalIndexLimits - returned array
!!
!! EXAMPLE
!!   Example 1. UG.
!!   For a 2d problem with a uniform grid block size of 8 
!!   With 4 blocks laid out in a square (2x2) grid
!! 
!!   globalIndexLimits(IAXIS) = 16
!!   globalIndexLimits(JAXIS) = 16
!!   globalIndexLimits(KAXIS) = 1  !because only 2d
!!
!!   Example 2. Paramesh
!!   For a 2d problem with block size of (8x8) and 
!!   maximum refinement level of 3, if the problem was
!!   initialized with one block then
!! 
!!   globalIndexLimits(IAXIS) = 32
!!   globalIndexLimits(JAXIS) = 32
!!   globalIndexLimits(KAXIS) = 1  !because only 2d
!!
!!   if problem was initialized with 2 blocks along IAXIS
!!   and 1 block along JAXIS then
!!
!!   globalIndexLimits(IAXIS) = 64
!!   globalIndexLimits(JAXIS) = 32
!!   globalIndexLimits(KAXIS) = 1  !because only 2d
!!
!!***

subroutine Grid_getGlobalIndexLimits(globalIndexLimits)

  implicit none
#include "constants.h"
  integer, dimension(MDIM), intent(out) :: globalIndexLimits
  globalIndexLimits=1
end subroutine Grid_getGlobalIndexLimits
