!!****f* source/Grid/Grid_receiveInputData
!!
!! NAME
!!
!!  Grid_receiveInputData
!!
!! SYNOPSIS
!!
!!  call Grid_receiveInputData(integer(IN) :: localNumBlocks,
!!                             integer(IN) :: alnblocks,
!!                             integer(IN) :: xx)
!!
!! DESCRIPTION 
!!
!!  Initializes grid arrays from arrays read by the I/O unit.
!!
!! ARGUMENTS  
!!
!!  localNumBlocks : the number of blocks on my processor.
!!
!!  alnblocks : the approximate number of local blocks on each
!!              processor if we give each processor an equal
!!              number of blocks.  Calculated from
!!              int(globalNumBlocks/meshNumProcs) + 1.
!!
!!  xx : an integer representing a cutoff point.  Processors
!!       less than this value are assigned alnblocks blocks and
!!       processors greater than or equal to this value are
!!       assigned lnblocks-1 blocks.
!!
!!***

subroutine Grid_receiveInputData(localNumBlocks, alnblocks, xx)
  implicit none
  integer, intent(IN) :: localNumBlocks, alnblocks, xx
end subroutine Grid_receiveInputData
