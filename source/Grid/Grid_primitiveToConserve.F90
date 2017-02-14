!!****f* source/Grid/Grid_primitiveToConserve
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

#include "constants.h"
#include "Flash.h"

subroutine Grid_primitiveToConserve(blkList,count,force)
  implicit none
  integer,intent(IN) :: count
  integer,dimension(count),intent(IN) :: blkList 
  logical,intent(IN) :: force
end subroutine Grid_primitiveToConserve
