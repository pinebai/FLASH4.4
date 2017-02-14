!!****if* source/Grid/GridMain/Chombo/AMR/Grid_restrictAllLevels
!!
!! NAME
!!  Grid_restrictAllLevels
!!
!! SYNOPSIS
!! 
!!  Grid_restrictAllLevels()
!!  
!! DESCRIPTION 
!!  Restricts the grid data to all refinement levels. Normally FLASH
!!  only evolves on the leaf blocks, calling this routine makes all
!!  levels have valid data.  This is mostly for visualization purposes to
!!  be able to look at different levels of resolution
!!  
!!  
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine Grid_restrictAllLevels()
  use Grid_data, ONLY : gr_convertToConsvdForMeshCalls
  use Grid_interface, ONLY : Grid_getListOfBlocks
  use Timers_interface, ONLY: Timers_start, Timers_stop
  use chombo_f_c_interface, ONLY : ch_restrict_all_levels
  implicit none
  integer, dimension(MAXBLOCKS) :: blkList
  integer :: count

  call Timers_start("restrictAll")

  if (gr_convertToConsvdForMeshCalls) then
     call Grid_getListOfBlocks(ALL_BLKS, blkList, count)
     call gr_primitiveToConserve(blkList, count)
  endif

  call ch_restrict_all_levels()

  if (gr_convertToConsvdForMeshCalls) then
     call gr_conserveToPrimitive(blkList, count, .FALSE.)
  endif

  call Timers_stop("restrictAll")
  
end subroutine Grid_restrictAllLevels
