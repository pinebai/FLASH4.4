!!****if* source/physics/materialProperties/Opacity/localAPI/op_LaguerreCoefficients
!!
!! NAME
!!
!!  op_LaguerreCoefficients
!!
!! SYNOPSIS
!!
!!  call op_LaguerreCoefficients (integer (in)  :: n,
!!                                real    (in)  :: alpha,
!!                                real    (out) :: LagA (1:n),
!!                                real    (out) :: LagB (1:n))
!!
!! DESCRIPTION
!!
!!  Generates the A and B coefficients defining the 3-term recursion relation for the
!!  generalized Laguerre polynomials:
!!
!!       alpha                 alpha            alpha
!!      Lag   (x) =  (x - A ) Lag   (x)  -  B  Lag   (x)
!!         i               i     i-1         i    i-2
!!
!!
!!       alpha                             alpha                               alpha
!!      Lag   (x) =  (x - 2i - alpha + 1) Lag   (x)  -  (i - 1)(i + alpha -1) Lag   (x)
!!         i                                 i-1                                 i-2
!!
!!  Thus we have:
!!
!!          A   = (2i + alpha - 1)              i = 1,2,...,n
!!           i
!!
!!          B   = (i - 1)(i + alpha -1)         i = 2,3,...,n
!!           i
!!
!! ARGUMENTS
!!
!!  n      : the maximum order of the polynomial wanted
!!  alpha  : the generalizing Laguerre factor
!!  LagA   : the polynomial A coefficients
!!  LagB   : the polynomial B coefficients
!!
!!***
subroutine op_LaguerreCoefficients (n, alpha, LagA, LagB)

  implicit none

  integer, intent (in)  :: n
  real,    intent (in)  :: alpha
  real,    intent (out) :: LagA (1:n)
  real,    intent (out) :: LagB (1:n)

  LagA = 0.0
  LagB = 0.0

  return
end subroutine op_LaguerreCoefficients
