!!****if* source/Grid/localAPI/gr_setBlockType
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

subroutine gr_setBlockType(blockID,type)
  implicit none

  integer, intent(IN) :: blockID, type

end subroutine gr_setBlockType
