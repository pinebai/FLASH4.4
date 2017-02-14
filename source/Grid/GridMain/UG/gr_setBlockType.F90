!!****if* source/Grid/GridMain/UG/gr_setBlockType
!!
!! NAME
!!  gr_setBlockType
!!
!! SYNOPSIS
!!
!!  gr_setBlockType(integer(IN)  :: blockID,
!!                     integer(IN)  :: type)
!!  
!! DESCRIPTION 
!!  
!!
!! ARGUMENTS 
!!
!!  blockID : the local blockid
!!  type    : the type to which block is to be set
!!
!!
!!
!! NOTES
!!
!!  gr_setBlockType is a mutator function that sets a block to the specified Type
!!  For example when tracing the path of a ray, a blocktype can be set to "TRAVERSE"
!!  if the ray passed through it.
!!
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine gr_setBlockType(blockID,type)
  use Grid_data, ONLY :  gr_blockType
  implicit none

  integer, intent(IN) :: blockID, type

  gr_blockType=type

end subroutine gr_setBlockType
