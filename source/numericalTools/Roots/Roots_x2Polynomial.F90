!!****f* source/numericalTools/Roots/Roots_x2Polynomial
!!
!! NAME
!!
!!  Roots_x2Polynomial
!!
!! SYNOPSIS
!!
!!  call Roots_x2Polynomial (real,    intent (in)  :: q1,
!!                           real,    intent (in)  :: q0,
!!                           integer, intent (out) :: nReal,
!!                           real,    intent (out) :: root (1:2,1:2))
!!
!! DESCRIPTION
!!
!!  Calculates all real + complex roots of the quadratic polynomial:
!!
!!                 x^2 + q1 * x + q0
!!
!!  The code checks internally, if rescaling of the coefficients is needed to
!!  avoid overflow.
!!
!!  The order of the roots is as follows:
!!
!!        1) For real roots, the order is according to their algebraic value
!!           on the number scale (largest positive first, largest negative last).
!!
!!        2) Since there can be only one complex conjugate pair root, no order
!!           is necessary.
!!
!! ARGUMENTS
!!
!!  q1         : coefficient of x term
!!  q0         : independent coefficient
!!  nReal      : number of real roots found
!!  root (n,1) : real part of n-th root
!!  root (n,2) : imaginary part of n-th root
!!
!! NOTES
!!
!!***

subroutine Roots_x2Polynomial (q1, q0,         &
                                        nReal, &
                                        root   )
  implicit none

  real   , intent (in)  :: q1, q0
  integer, intent (out) :: nReal
  real   , intent (out) :: root (1:2,1:2)

  nReal = 0
  root (1:2,1:2) = 0.0

  return
end subroutine Roots_x2Polynomial
