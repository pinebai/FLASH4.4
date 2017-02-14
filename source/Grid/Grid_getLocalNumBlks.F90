!!****f* source/Grid/Grid_getLocalNumBlks
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
  implicit none
  integer,intent(out) :: numBlocks

  !assign default value for stub
  numBlocks = 0
  return
end subroutine Grid_getLocalNumBlks
