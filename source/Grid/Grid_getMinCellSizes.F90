!!****f* source/Grid/Grid_getMinCellSizes
!!
!! NAME
!!  Grid_getMinCellSizes
!!
!! SYNOPSIS
!!
!!  Grid_getMinCellSizes(real (OUT)  :: minCellSizes(MDIM))
!!               
!!  
!! DESCRIPTION 
!!
!!  Returns the smallest possible cell size in a simulation for all dimensions.
!!
!!
!! ARGUMENTS
!!
!!  minCellSizes - returned array
!!
!!***

#include "constants.h"

subroutine Grid_getMinCellSizes(minCellSizes)
  implicit none
  real, dimension(MDIM), intent(OUT) :: minCellSizes
  minCellSizes(1:MDIM) = 0
end subroutine Grid_getMinCellSizes
