!!****if* source/Grid/GridMain/paramesh/Grid_getGlobalIndexLimits
!!
!! NAME
!!  Grid_getGlobalIndexLimits
!!
!! SYNOPSIS
!!
!!  Grid_getGlobalIndexLimits(integer :: globalIndexLimits(MDIM))
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
!!  In an adaptive mesh, the highest index is returned as if
!!  the entire mesh was fully refined.
!!
!! ARGUMENTS
!!  globalIndexLimits - returned array
!!
!! EXAMPLE
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

  use Grid_data, ONLY : gr_nblockX, gr_nblockY, gr_nblockZ
  use tree, ONLY : lrefine_max

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, dimension(MDIM), intent(OUT) :: globalIndexLimits
  integer, dimension(MDIM) :: maxBlocksAlongDimension

  maxBlocksAlongDimension(1:MDIM) = 1   !Required when MDIM > NDIM.
  maxBlocksAlongDimension(1:NDIM) = 2 ** (lrefine_max-1)

  globalIndexLimits(IAXIS) = NXB * gr_nblockX * maxBlocksAlongDimension(IAXIS)
  globalIndexLimits(JAXIS) = NYB * gr_nblockY * maxBlocksAlongDimension(JAXIS)
  globalIndexLimits(KAXIS) = NZB * gr_nblockZ * maxBlocksAlongDimension(KAXIS)

end subroutine Grid_getGlobalIndexLimits
