!!****if* source/Grid/GridMain/paramesh/Grid_getBlkRefineLevel
!!
!!
!! NAME
!!  Grid_getBlkRefineLevel
!!
!! SYNOPSIS
!!
!!
!!  Grid_getBlkRefineLevel(integer(IN)  :: blockID,
!!                      integer(OUT) :: refineLevel)
!!  
!! DESCRIPTION 
!!  Get the refinement level of a given block as denoted by blockId
!!
!! ARGUMENTS
!!  blockID - the local block number
!!  refineLevel - returned value, refinement level of block
!!
!!***

subroutine Grid_getBlkRefineLevel(blockID,refineLevel)

  use tree, ONLY : lrefine

  implicit none
  integer,intent(in) :: blockID
  integer,intent(out) :: refineLevel

  refineLevel = lrefine(blockId)

  return
end subroutine Grid_getBlkRefineLevel














