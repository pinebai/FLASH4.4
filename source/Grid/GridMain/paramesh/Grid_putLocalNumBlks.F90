!!****if* source/Grid/GridMain/paramesh/Grid_putLocalNumBlks
!!
!! NAME
!!  Grid_putLocalNumBlks
!!
!! SYNOPSIS
!!
!!  Grid_putLocalNumBlks(integer(IN) :: numBlocks)
!!  
!! DESCRIPTION 
!!  Put the number of local blocks on a processor 
!!  Only used in restart capabilities
!!
!! ARGUMENTS 
!!  numBlocks : the local number of blocks currently in use,
!!              supplied by the caller 
!! 
!! 
!!
!!***


subroutine Grid_putLocalNumBlks(numBlocks)

  use tree, ONLY : lnblocks

  implicit none
  integer,intent(in) :: numBlocks

  lnblocks = numBlocks

  return
end subroutine Grid_putLocalNumBlks
