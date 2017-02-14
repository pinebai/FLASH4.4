!!****if* source/Grid/GridMain/paramesh/Grid_markBlkRefine
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
!! NOTES
!!
!! Every unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. For Grid unit these variables begin with "gr_"
!! like, gr_meshMe or gr_eosMode, and are stored in fortran
!! module Grid_data (in file Grid_data.F90. The other variables
!! are local to the specific routines and do not have the prefix "gr_"
!!
!! 
!!
!!***

subroutine Grid_markBlkRefine(block,mark)
  use tree, ONLY : refine
implicit none
  integer, intent(IN) :: block
  logical, intent(IN) :: mark

  refine(block)   = mark

end subroutine Grid_markBlkRefine
