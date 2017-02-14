!!****if* source/Grid/GridMain/Grid_getMinCellSizes
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
  use Grid_data, ONLY : gr_minCellSizes
  implicit none
  real, dimension(MDIM), intent(OUT) :: minCellSizes
  minCellSizes(1:MDIM) = gr_minCellSizes(1:MDIM)
end subroutine Grid_getMinCellSizes
