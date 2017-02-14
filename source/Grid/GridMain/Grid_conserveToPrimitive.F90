!!****if* source/Grid/GridMain/Grid_conserveToPrimitive
!!
!! NAME
!!
!!  Grid_conserveToPrimitive
!!
!!
!! SYNOPSIS
!!
!!  Grid_conserveToPrimitive(integer(in) :: blkList(count),
!!                           integer(in) :: count,
!!                           logical(in) :: allCells,
!!                           logical(in) :: force)
!!
!!
!! DESCRIPTION
!!
!!  Calls gr_conserveToPrimitive
!!
!!
!! ARGUMENTS
!! 
!!   blkList - integer list of blocks to be operated on
!!
!!   count - number of blocks in the blkList
!!
!!   allCells - act on all cells, including guardcells, if .TRUE.,
!!              otherwise only modify interior cells.
!!
!!   force - whether to force conversion
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine Grid_conserveToPrimitive(blkList,count,allCells,force)

#ifndef FLASH_GRID_UG
  use Grid_data, ONLY: gr_convertToConsvdForMeshCalls
#endif

  implicit none
  integer,intent(IN) :: count
  integer,dimension(count),intent(IN) :: blkList 
  logical,intent(IN) :: allCells, force
  logical :: tempSwap

#ifndef FLASH_GRID_UG
  tempSwap = (force .and. (.not.gr_convertToConsvdForMeshCalls))
  if (tempSwap) then
     gr_convertToConsvdForMeshCalls = .true.
  end if

  call gr_conserveToPrimitive(blkList,count,allCells)

  if (tempSwap) then
     gr_convertToConsvdForMeshCalls = .false.
  end if
#endif

end subroutine Grid_conserveToPrimitive
