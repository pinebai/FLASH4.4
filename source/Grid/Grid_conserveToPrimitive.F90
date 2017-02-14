!!****f* source/Grid/Grid_conserveToPrimitive
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
  implicit none
  integer,intent(IN) :: count
  integer,dimension(count),intent(IN) :: blkList 
  logical,intent(IN) :: allCells, force
end subroutine Grid_conserveToPrimitive
