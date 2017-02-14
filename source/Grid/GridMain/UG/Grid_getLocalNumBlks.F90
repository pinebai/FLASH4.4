!!****if* source/Grid/GridMain/UG/Grid_getLocalNumBlks
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

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif


subroutine Grid_getLocalNumBlks(numBlocks)
implicit none
  integer,intent(out) :: numBlocks
  numBlocks = 1
  return
end subroutine Grid_getLocalNumBlks
