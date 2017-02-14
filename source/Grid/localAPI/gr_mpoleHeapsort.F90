!!****if* source/Grid/localAPI/gr_mpoleHeapsort
!!
!! NAME
!!
!!  gr_mpoleHeapsort
!!
!! SYNOPSIS
!!
!!  gr_mpoleHeapsort (integer, intent(in)    :: nElements,
!!                    real,    intent(inout) :: Vector )
!!
!! DESCRIPTION
!!
!!  Applies the heapsort algorithm on a Vector of size 'nElements'.
!!  When exiting, the vector elements are ordered in increasing order:
!!
!!           Vector (i+1) >= Vector (i)      1 <= i <= nElements
!!
!!  The algorithm scales as O(N*log2(N)) for the worst case scenario.
!!
!! ARGUMENTS
!!
!!  nElements : number of vector elements
!!  Vector    : the vector itself (in: not ordered, out: ordered)
!!
!!***

subroutine gr_mpoleHeapsort (nElements,Vector)

  implicit none
  
  integer, intent (in)    :: nElements
  real,    intent (inout) :: Vector (1:nElements)

  Vector = 0.0

  return
end subroutine gr_mpoleHeapsort
