!!****f* source/Grid/Grid_getBlkRefineLevel
!!
!! NAME
!!  Grid_getBlkRefineLevel
!!
!! SYNOPSIS
!!
!!
!!  Grid_getBlkRefineLevel(integer(IN)  :: blockID,
!!                         integer(OUT) :: refineLevel)
!!  
!! DESCRIPTION 
!!  Get the refinement level of a given block as denoted by blockID.
!!
!! ARGUMENTS
!!  blockID - the local block number
!!  refineLevel - returned value, refinement level of block
!!
!!***

subroutine Grid_getBlkRefineLevel(blockID, refineLevel)
  implicit none
  integer,intent(in) :: blockID
  integer,intent(out) :: refineLevel
  refineLevel=1
  return
end subroutine Grid_getBlkRefineLevel














