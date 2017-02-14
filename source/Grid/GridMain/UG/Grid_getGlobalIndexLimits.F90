!!****if* source/Grid/GridMain/UG/Grid_getGlobalIndexLimits
!!
!! NAME
!!  Grid_getGlobalIndexLimits
!!
!! SYNOPSIS
!!
!!  Grid_getGlobalIndexLimits(integer(OUT) :: globalIndexLimits(MDIM))
!!  
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
!! ARGUMENTS
!!  globalIndexLimits - returned array
!!
!! EXAMPLE
!!   For a 2d problem with a uniform grid block size of 8 
!!   With 4 blocks laid out in a square (2x2) grid
!! 
!!   globalIndexLimits(IAXIS) = 16
!!   globalIndexLimits(JAXIS) = 16
!!   globalIndexLimits(KAXIS) = 1  !because only 2d
!!
!!
!!***



subroutine Grid_getGlobalIndexLimits(globalIndexLimits)

  use Grid_data, ONLY : gr_gIndexSize

  implicit none

#include "constants.h"

  integer, dimension(MDIM), intent(OUT) :: globalIndexLimits

  globalIndexLimits = gr_gIndexSize

end subroutine Grid_getGlobalIndexLimits
