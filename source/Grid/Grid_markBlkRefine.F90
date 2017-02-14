!!****f* source/Grid/Grid_markBlkRefine
!!
!! NAME
!!  Grid_markBlkRefine
!!
!! SYNOPSIS
!!
!!  Grid_markBlkRefine(integer(IN) :: block,
!!                       logical(IN) :: mark)
!!  
!! DESCRIPTION 
!!  Mark a specific block for refinement.
!!  The paramesh version of Grid implements this function, whereas
!!  UG uses the stub.
!!
!! ARGUMENTS
!!  block : local Block id
!!  mark : mark the block if it should be refined
!! 
!!
!! 
!!
!!***

subroutine Grid_markBlkRefine(block,mark)
implicit none
  integer, intent(IN) :: block
  logical, intent(IN) :: mark

end subroutine Grid_markBlkRefine
