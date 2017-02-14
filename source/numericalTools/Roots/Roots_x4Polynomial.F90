!!****f* source/numericalTools/Roots/Roots_x4Polynomial
!!
!! NAME
!!
!!  Roots_x4Polynomial
!!
!! SYNOPSIS
!!
!!  call Roots_x4Polynomial (real,              intent (in)  :: q3,
!!                           real,              intent (in)  :: q2,
!!                           real,              intent (in)  :: q1,
!!                           real,              intent (in)  :: q0,
!!                           integer,           intent (out) :: nReal,
!!                           real,              intent (out) :: root (1:4,1:2),
!!                           logical, optional, intent (in)  :: printInfo,
!!                           integer, optional, intent (in)  :: printUnit)
!!
!! DESCRIPTION
!!
!!  Calculates all real + complex roots of the quartic polynomial:
!!
!!                 x^4 + q3 * x^3 + q2 * x^2 + q1 * x + q0
!!
!!  An option for printing a detailed info about the intermediate stages in solving
!!  the quartic is available. Since the code has not yet been extensively tested,
!!  this enables a detailed check in case something went wrong and the roots obtained
!!  are not proper.
!!
!!  The quartic root solver can handle any size of quartic coefficients and there is
!!  no danger of overflow, due to proper rescaling of the quartic polynomial.
!!
!!  The order of the roots is as follows:
!!
!!        1) For real roots, the order is according to their algebraic value
!!           on the number scale (largest positive first, largest negative last).
!!
!!        2) For complex conjugate pair roots, the order is according to the
!!           algebraic value of their real parts (largest positive first). If
!!           the real parts are equal, the order is according to the algebraic
!!           value of their imaginary parts (largest first).
!!
!!        3) All real roots preceede the complex ones.
!!
!! ARGUMENTS
!!
!!  q3         : coefficient of x^3 term
!!  q2         : coefficient of x^2 term
!!  q1         : coefficient of x term
!!  q0         : independent coefficient
!!  nReal      : number of different real roots found
!!  root (n,1) : real part of n-th root
!!  root (n,2) : imaginary part of n-th root
!!  printInfo  : if given and true, detailed info will be printed about intermediate stages
!!  printUnit  : the unit ID, where the info will be printed
!!
!! NOTES
!!
!!  If omitting 'printUnit' but specifying 'printInfo' as true, the info will be printed
!!  on default output monitor. Giving just 'printUnit' does not do anything, unless 'printInfo'
!!  is given and set to true.
!!
!!***

subroutine Roots_x4Polynomial (q3, q2, q1, q0,        &
                                               nReal, &
                                               root,  &
                               printInfo,             &
                               printUnit              )
  implicit none

  real   ,           intent (in)  :: q3, q2, q1, q0
  integer,           intent (out) :: nReal
  real   ,           intent (out) :: root (1:4,1:2)
  logical, optional, intent (in)  :: printInfo
  integer, optional, intent (in)  :: printUnit

  nReal = 0
  root (1:4,1:2) = 0.0

  return
end subroutine Roots_x4Polynomial
