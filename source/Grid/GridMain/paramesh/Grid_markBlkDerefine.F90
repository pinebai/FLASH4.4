!!****if* source/Grid/GridMain/paramesh/Grid_markBlkDerefine
!!
!! NAME
!!  Grid_markBlkDerefine
!!
!! SYNOPSIS
!!
!!
!! SYNOPSIS
!!
!!  Grid_markBlkDerefine(integer(IN) :: blockID,
!!                       logical(IN) :: mark)
!!  
!! DESCRIPTION 
!!  Mark a specific block for derefinement.
!!  The paramesh version of Grid implements this function, whereas
!!  UG uses the stub.
!!
!! ARGUMENTS
!!  blockID : local Block id
!!  mark : mark the block if it should be derefined
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
!!*** 
subroutine Grid_markBlkDerefine(blockID,mark)
  use tree, ONLY : derefine
  implicit none
  integer, intent(IN) :: blockID
  logical, intent(IN) :: mark

  derefine(blockID)   = mark

end subroutine Grid_markBlkDerefine
