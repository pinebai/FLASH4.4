!!****f* source/Grid/Grid_getBlkType
!!
!! NAME
!!  Grid_getBlkType
!!
!! SYNOPSIS
!!
!!
!!  Grid_getBlkType(integer(IN)  :: blockID,
!!                  integer(OUT) :: blkType)
!!  
!! DESCRIPTION 
!!  Get the type of block. This is relevant only in AMR mode.
!!  Valid types are described in constants.h, and include
!!  LEAF, PARENT and ANCESTOR for paramesh.
!!
!! ARGUMENTS
!!  blockID - the local block number
!!  blkType - returned value
!!
!!***

subroutine Grid_getBlkType(blockID, blkType)
  implicit none
#include "constants.h"
  integer,intent(in) :: blockID
  integer,intent(out) :: blkType
  blkType=LEAF
  return
end subroutine Grid_getBlkType














