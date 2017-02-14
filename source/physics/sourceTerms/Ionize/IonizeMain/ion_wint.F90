!!****if* source/physics/sourceTerms/Ionize/IonizeMain/ion_wint
!!
!! NAME
!!  
!!  ion_wint
!!
!! SYNOPSIS
!! 
!!  ion_wint(real(IN)     :: xaux(n),
!!           integer(IN)  :: n,
!!           real(IN)     :: x,
!!           integer(OUT) :: j)
!!
!!  
!! DESCRIPTION
!!
!! Routine that gives the index of xaux immediately before of x
!!
!! ARGUMENTS
!!
!! n:    size of xaux.
!! xaux: array of possible values.
!! x:    number to look for.
!! j:    index you want.
!!
!!***

subroutine ion_wint(xaux,n,x,j)
  implicit none

  integer, intent(IN) :: n
  real, dimension(n), intent(IN) :: xaux
  integer, intent(OUT) :: j
  real, intent(IN) :: x

  integer :: jl
  integer ::  ju
  integer :: jm 
  logical :: A,B

  jl = 0 
  ju = n+1
  do while (ju-jl .gt. 1)   
     jm = (ju+jl )/ 2
     A = ((xaux(n) - xaux(1)) > 0.0)
     B = (x - xaux(jm) > 0.0)
     if (A .eqv. B) then
        jl = jm
     else
        ju = jm
     end if
  end do
  j = jl
  j = min(j, n-1)
  j = max(j, 1)

  return
end subroutine ion_wint
