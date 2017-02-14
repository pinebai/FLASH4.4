!!****if* source/flashUtilities/interpolation/oneDim/ut_polint
!!
!!  NAME
!!    ut_polint
!!
!!  SYNOPSIS 
!!    call ut_polint( real(in)    :: xa(*),
!!                    real(in)    :: ya(*),
!!                    real(in)    :: n,
!!                    real(in)    :: x,
!!                    real(out)   :: y,
!!                    real(inout) :: dy)
!!
!!  DESCRIPTION
!!    Polynomial interpolation routine with the same interface
!!    as POLINT from Numerical Recipies.
!!
!!    Only linear and quadratic interpolation are actually implemented.
!!
!!  FUNCTION
!!    Given n-element arrays (xa,ya) of locations and values, and an
!!    interpolation point x, returns a polynomial interpolant
!!    (or extrapolant -- no range checking done!) of the function in y.
!!
!!  ARGUMENTS
!!
!!      xa - array of real numbers
!!      ya - array of real numbers
!!      n - length of arrays. Number of input location/value pairs used
!!          as input by the interpolation.
!!          This determines the order of interpolation.
!!      x - value to be interpolated
!!      y - interpolated value
!!      dy - error estimate of y.
!!           Provided for call compatibility with POLINT, left unchanged
!!           by this implementation.
!!
!!  RESULT
!!      returns result in y
!!
!!  HISTORY
!!
!!    Created as a wrapper            - KW 2007-06-07
!!***
subroutine ut_polint(xa, ya, n, x, y, dy)

  use Driver_interface, ONLY : Driver_abortFlash
  use ut_interpolationInterface, ONLY : ut_quadraticInterpol

  implicit none
  real, intent(in)    :: x           ! where function is to be eval'd
  integer,intent(in)  :: n           ! number of input location/value pairs
  real, intent(in)    :: ya(*), xa(*)     ! function values and locations
  real, intent(out)   :: y              ! function evaluated at x
  real, intent(inout) :: dy             ! ignored and unchanged


  select case (n)
  case(3)
     call ut_quadraticInterpol(x=xa, y=ya, xint=x, yint=y)
  case(2)
     if (abs(x-xa(1)) .LE. abs(x-xa(2))) then
        y = ya(1)  +  (ya(2) - ya(1)) / (xa(2) - xa(1))  *  (x - xa(1))
     else
        y = ya(2)  +  (ya(2) - ya(1)) / (xa(2) - xa(1))  *  (x - xa(2))
     end if
  case default
     call Driver_abortFlash('ut_polint only implements linear and parabolic interpolation!')
  end select

  return
end subroutine ut_polint
