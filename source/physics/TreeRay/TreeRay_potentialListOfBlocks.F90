!!****f* source/physics/TreeRay/TreeRay_potentialListOfBlocks
!!
!!  NAME 
!!
!!     TreeRay_potentialListOfBlocks
!!
!!  SYNOPSIS
!!
!!  TreeRay_potentialListOfBlocks(integer(IN) :: blockCount,
!!                                integer(IN) :: blockList(blockCount))
!!
!!  DESCRIPTION 
!!
!! ARGUMENTS
!!
!!   blockCount   : The number of blocks in the list
!!   blockList(:) : The list of blocks on which to calculate potential
!!
!! SIDE EFFECTS
!!
!!
!!
!!***


subroutine TreeRay_potentialListOfBlocks(blockCount,blockList)


  implicit none

  integer,intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN) :: blockList

!=========================================================================

  
  return
end subroutine TreeRay_potentialListOfBlocks

