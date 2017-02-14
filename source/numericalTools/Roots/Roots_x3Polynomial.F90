!!****f* source/numericalTools/Roots/Roots_x3Polynomial
!!
!! NAME
!!
!!  Roots_x3Polynomial
!!
!! SYNOPSIS
!!
!!  call Roots_x3Polynomial (real,              intent (in)  :: c2,
!!                           real,              intent (in)  :: c1,
!!                           real,              intent (in)  :: c0,
!!                           integer,           intent (out) :: nReal,
!!                           real,              intent (out) :: root (1:3,1:2),
!!                           logical, optional, intent (in)  :: printInfo,
!!                           integer, optional, intent (in)  :: printUnit)
!!
!! DESCRIPTION
!!
!!  Calculates all real + complex roots of the cubic polynomial:
!!
!!                 x^3 + c2 * x^2 + c1 * x + c0
!!
!!  The first real root (which always exists) is obtained using an optimized
!!  Newton-Raphson scheme. The other remaining roots are obtained through
!!  composite deflation to a quadratic.
!!
!!  The cubic root solver can handle any size of cubic coefficients and there is
!!  no danger of overflow due to proper rescaling of the cubic polynomial.
!!
!!  The order of the roots is as follows:
!!
!!        1) For real roots, the order is according to their algebraic value
!!           on the number scale (largest positive first, largest negative last).
!!
!!        2) Since there can be only one complex conjugate pair root, no order
!!           is necessary.
!!
!!        3) All real roots preceede the complex ones.
!!
!! ARGUMENTS
!!
!!  c2         : coefficient of x^2 term
!!  c1         : coefficient of x term
!!  c0         : independent coefficient
!!  nReal      : number of different real roots found
!!  root (n,1) : real part of n-th root
!!  root (n,2) : imaginary part of n-th root
!!  printInfo  : if given and true, detailed info will be printed about intermediate stages
!!  printUnit  : the unit ID, where the info will be printed
!!
!! NOTES
!!
!!  Only passing both printInfo AND printUnit will result in printing out the info.
!!  Giving only one of them results in no printing action.
!!
!!***

subroutine Roots_x3Polynomial (c2, c1, c0,         &
                                            nReal, &
                                            root,  &
                               printInfo,          &
                               printUnit           )
  implicit none

  real   ,           intent (in)  :: c2, c1, c0
  integer,           intent (out) :: nReal
  real   ,           intent (out) :: root (1:3,1:2)
  logical, optional, intent (in)  :: printInfo
  integer, optional, intent (in)  :: printUnit

  nReal = 0
  root (1:3,1:2) = 0.0

  return
end subroutine Roots_x3Polynomial
