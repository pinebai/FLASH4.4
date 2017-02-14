!!****if* source/Grid/GridMain/Chombo/Grid_getLocalNumBlks
!!
!! NAME
!!  Grid_getLocalNumBlks
!!
!! SYNOPSIS
!!
!!  Grid_getLocalNumBlks(numBlocks)
!!  Grid_getLocalNumBlks(int(out))
!!  
!! DESCRIPTION 
!!  Get the number of local blocks on a processor 
!!
!!
!! CREATED : 05/18/04, by AD
!! 
!! 
!! 
!! 
!!
!!***

subroutine Grid_getLocalNumBlks(numBlocks)
  use chombo_f_c_interface, ONLY : ch_get_num_blocks
  implicit none
  integer,intent(out) :: numBlocks

  call ch_get_num_blocks(numBlocks)
end subroutine Grid_getLocalNumBlks
