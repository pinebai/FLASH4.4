!!****if* source/Grid/GridMain/Grid_primitiveToConserve
!!
!! NAME
!!
!!  Grid_primitiveToConserve
!!
!!
!! SYNOPSIS
!!
!!  Grid_primitiveToConserve(integer(in) :: blkList(count),
!!                           integer(in) :: count,
!!                           logical(in) :: force)
!!
!!
!! DESCRIPTION
!!
!!  Calls gr_primitiveToConserve
!!
!!
!! ARGUMENTS
!! 
!!   blkList - integer list of blocks to be operated on
!!
!!   count - number of blocks in the blkList
!!
!!   force - whether to force conversion
!!
!!***

#include "Flash.h"

subroutine Grid_primitiveToConserve(blkList,count,force)

#ifndef FLASH_GRID_UG
  use Grid_data, ONLY: gr_convertToConsvdForMeshCalls
#endif

  implicit none
  integer,intent(IN) :: count
  integer,dimension(count),intent(IN) :: blkList 
  logical,intent(IN) :: force
  logical :: tempSwap

#ifndef FLASH_GRID_UG
  tempSwap = (force .and. (.not.gr_convertToConsvdForMeshCalls))
  if (tempSwap) then
     gr_convertToConsvdForMeshCalls = .true.
  end if

  call gr_primitiveToConserve(blkList,count)

  if (tempSwap) then
     gr_convertToConsvdForMeshCalls = .false.
  end if
#endif

end subroutine Grid_primitiveToConserve
