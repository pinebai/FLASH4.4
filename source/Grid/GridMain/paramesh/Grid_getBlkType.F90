!!****if* source/Grid/GridMain/paramesh/Grid_getBlkType
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
!!  Get the type of block. 
!!  Valid types are defined in constants.h, and include
!!  LEAF, PARENT and ANCESTOR for paramesh.
!!
!! ARGUMENTS
!!  blockID - the local block number
!!  blkType - returned value
!!
!!***

subroutine Grid_getBlkType(blockID, blkType)
  use tree, ONLY : nodetype
  implicit none
#include "constants.h"
  integer,intent(in) :: blockID
  integer,intent(out) :: blkType
  blkType=nodetype(blockID)
  return
end subroutine Grid_getBlkType














