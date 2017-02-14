!!****if* source/Grid/GridMain/UG/Grid_getBlkRefineLevel
!!
!! NAME
!!  Grid_getBlkRefineLevel
!!
!! SYNOPSIS
!!
!!  Grid_getBlkRefineLevel(integer(IN)  :: blockID,
!!                         integer(OUT) :: refineLevel)
!!  
!! DESCRIPTION 
!!  Get the refinement level of a given block as denoted by blockID.
!!  For the UG this is always 1.
!!
!! ARGUMENTS
!!  blockID - the local block number
!!  refineLevel - returned value, refinement level of block
!!
!! 
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine Grid_getBlkRefineLevel(blockID,refineLevel)
implicit none
  integer,intent(in) :: blockID
  integer,intent(out) :: refineLevel

  refineLevel = 1

  return
end subroutine Grid_getBlkRefineLevel














