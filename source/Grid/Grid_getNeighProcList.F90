!!****f* source/Grid/Grid_getNeighProcList
!!
!! NAME
!!  Grid_getNeighProcList
!!
!! SYNOPSIS
!!
!!  call Grid_getNeighProcList(logical(IN)             :: includeMyProc, 
!!                             integer(INOUT), pointer :: neighProcList(:),
!!                             integer(OUT)            :: numNeigh)
!!
!! DESCRIPTION 
!!
!!  Creates a pointer array containing the neighboring processor IDs of all
!!  LEAF blocks on this processor.
!!
!! ARGUMENTS 
!!  
!!
!!  includeMyProc - Whether the array should include my processor ID.
!!
!!  neighProcList - The processor IDs of all neighboring LEAF blocks. 
!!
!!  numNeigh      - The number of entries in neighProcList.
!!
!!
!! NOTES
!!
!!  Currently only implemented for Paramesh4 Grid implementations.
!!
!!  It is the users resposibility to
!!    1. deallocate neighProcList
!!    2. obtain a new neighProcList when the mesh changes
!!
!!***

subroutine Grid_getNeighProcList(includeMyProc, neighProcList, numNeigh)
  implicit none
  logical, intent(IN) :: includeMyProc
  integer, dimension(:), pointer :: neighProcList
  integer, intent(OUT) :: numNeigh
  numNeigh = 0
end subroutine Grid_getNeighProcList
