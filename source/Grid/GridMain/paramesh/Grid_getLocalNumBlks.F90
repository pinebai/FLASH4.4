!!****if* source/Grid/GridMain/paramesh/Grid_getLocalNumBlks
!!
!! NAME
!!  Grid_getLocalNumBlks
!!
!! SYNOPSIS
!!
!!  Grid_getLocalNumBlks(integer(OUT) :: numBlocks)
!!  
!! DESCRIPTION 
!!  Get the number of local blocks on a processor 
!!
!! ARGUMENTS
!!  numBlocks : The number of blocks currently in use on myProcessor
!!
!!***


subroutine Grid_getLocalNumBlks(numBlocks)
  
  use tree, ONLY : lnblocks
  
  implicit none
  
  integer,intent(out) :: numBlocks
  
  numBlocks = lnblocks
  
  return
end subroutine Grid_getLocalNumBlks
