!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpoleHeapsort
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

  integer :: halfArray
  integer :: i,j,n

  real    :: Vectorn
!
!
!    ...The first step of the heapsort: arrange the initial Vector array
!       into heap form.
!
!
  halfArray = nElements / 2

  do n = halfArray,1,-1
     Vectorn = Vector (n)
     i = n
     j = n + n
     do while (j <= nElements)
        if (j < nElements) then
            if (Vector (j) < Vector (j+1)) then
                j = j + 1
            end if
        end if
        if (Vectorn < Vector (j)) then
            Vector (i) = Vector (j)
            i = j
            j = j + j
        else
            j = nElements + 1
        end if
     end do
     Vector (i) = Vectorn
  end do
!
!
!    ...The second step of the heapsort: promotion and re-heapify.
!
!
  do n = nElements,2,-1
     Vectorn = Vector (n)
     Vector (n) = Vector (1)
     i = 1
     j = 2
     do while (j <= n-1)
        if (j < n-1) then
            if (Vector (j) < Vector (j+1)) then
                j = j + 1
            end if
        end if
        if (Vectorn < Vector (j)) then
            Vector (i) = Vector (j)
            i = j
            j = j + j
        else
            j = n
        end if
     end do
     Vector (i) = Vectorn
  end do
!
!
!    ...Ready!
!
!
  return
end subroutine gr_mpoleHeapsort
