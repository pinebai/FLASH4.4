!!****f* source/Grid/Grid_markBlkDerefine
!!
!! NAME
!!  Grid_markBlkDerefine
!!
!! SYNOPSIS
!!
!!  Grid_markBlkDerefine(integer(IN) :: block,
!!                       logical(IN) :: mark)
!!  
!! DESCRIPTION 
!!  Mark a specific block for derefinement.
!!  The paramesh version of Grid implements this function, whereas
!!  UG uses the stub.
!!
!! ARGUMENTS
!!  block : local Block id
!!  mark : mark the block if it should be derefined
!! 
!! 
!!
!!***

subroutine Grid_markBlkDeRefine(block,mark)
implicit none
  integer, intent(IN) :: block
  logical, intent(IN) :: mark

end subroutine Grid_markBlkDeRefine
