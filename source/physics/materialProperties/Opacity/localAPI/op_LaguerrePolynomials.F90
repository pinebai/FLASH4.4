!!****if* source/physics/materialProperties/Opacity/localAPI/op_LaguerrePolynomials
!!
!! NAME
!!
!!  op_LaguerrePolynomials
!!
!! SYNOPSIS
!!
!!  call op_LaguerrePolynomials (integer (in)  :: n,
!!                               integer (in)  :: nLag,
!!                               real    (in)  :: alpha,
!!                               real    (in)  :: X,
!!                               real    (out) :: Lag (1:nLag))
!!
!! DESCRIPTION
!!
!!  Generates up to n-th order rescaled generalized Laguerre polynomials, defined
!!  by the following recursion:
!!
!!       alpha        X - 2i - alpha + 1   alpha           (i - 1)(i + alpha -1)      alpha
!!      Lag   (X) =  -------------------- Lag   (X)  -  ---------------------------- Lag   (X)
!!         i              X +/- 2i           i-1         (X +/- 2i)(X +/- 2[i - 1])     i-2
!!
!!  where
!!
!!             +/-  means:  take - sign for X < 1.0 and + sign otherwise
!!
!!
!!  The recursion is extremely stable in the forward direction, retaining full precision as
!!  the i-index increases. The 0-th order polynomial value of 1. is not returned or stored!
!!
!!
!!  Motivation for rescaling the polynomials
!!  ----------------------------------------
!!
!!  Consider the standard generalized Laguerre polynomial recursion relation:
!!
!!
!!          alpha                            alpha                               alpha
!!         L     (X) = (X - 2i - alpha + 1) L     (X)  -  (i - 1)(i + alpha -1) L     (X)
!!          i                                i-1                                 i-2
!!
!!
!!                                        alpha
!!  If (X >> i) the polynomials scale as L     (X) ->  R^i
!!                                        i
!!
!!                                        alpha
!!  If (X << i) the polynomials scale as L     (X) ->  2^i * i!
!!                                        i
!!
!!  The introduction of the scaling factor:
!!
!!                          k=i
!!                 P   =   Prod  (X +/- 2k)
!!                  i       k=1
!!
!!  renders the recursion much more stable against computational under- and overflow.
!!  Note that the ideal rescaling factor involving only (X - 2k) factors would do no
!!  good due to zeroes at arguments for which X = 2k and very small numbers in the
!!  neighborhood of X = 2k.
!!
!!  When using the rescaled generalized Laguerre polynomials to establish moments
!!  of the standard generalized Laguerre polynomials over a weight function, the scaling
!!  factor must be considered as well!
!!
!! ARGUMENTS
!!
!!  n     : the maximum order of the polynomial wanted
!!  nLag  : the declared dimension for the array holding the Laguerre polynomials
!!  alpha : the generalizing factor
!!  X     : the polynomial argument
!!  Lag   : the polynomials for i = 1,...,n
!!
!!***
subroutine op_LaguerrePolynomials (n, nLag, alpha, X, Lag)

  implicit none

  integer, intent (in)  :: n
  integer, intent (in)  :: nLag
  real,    intent (in)  :: alpha
  real,    intent (in)  :: X
  real,    intent (out) :: Lag (1:nLag)

  Lag = 0.0

  return
end subroutine op_LaguerrePolynomials
