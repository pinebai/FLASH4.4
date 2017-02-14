!!****if* source/physics/materialProperties/Opacity/localAPI/op_LaguerreQuadratureRule
!!
!! NAME
!!
!!  op_LaguerreQuadratureRule
!!
!! SYNOPSIS
!!
!!  call op_LaguerreQuadratureRule (integer (in)  :: nRoots,
!!                                  real    (in)  :: beta,
!!                                  real    (in)  :: T,
!!                                  real    (out) :: Roots,
!!                                  real    (out) :: Weights)
!!
!! DESCRIPTION
!!
!!  Generates 'nRoots' number of roots and weights for integration of Laguerre-type integrands:
!!
!!                              T
!!                            /
!!                           |         beta  -x
!!             Integral =    |  F (x) x     e   dx
!!                           |
!!                          /
!!                         0
!!
!!  using alpha-generalized Laguerre polynomials as auxilliary polynomials. The integral value
!!  will be exact, if F(x) is a polynomial of order (2*nRoots - 1) and very well approximated,
!!  if F(x) can be represented with high accuracy as a polynomial of such order.
!!
!!  Main sequence of steps
!!  ----------------------
!!
!!     1) Establish the generalized Laguerre polynomials for the integration limits T and 0
!!     2) Calculate the associated modified moments over the Laguerre-type weight x^beta e^(-x)
!!     3) Normalize the modified moments (0-th moment = 1)
!!     4) Set up terminal matrix (Sack & Donovan and Wheeler algorithm)
!!     5) Find roots and weights from the terminal matrix using the Golub-Welsch algorithm
!!
!!
!! ARGUMENTS
!!
!!  nRoots  : the number of roots and weights wanted
!!  beta    : the exponential for the Laguerre-type weight function (whole number for now!)
!!  T       : the upper integration limit (must be > 0)
!!  Roots   : the roots for i = 1,2,...,nRoots
!!  Weights : the weights for i = 1,2,...,nRoots
!!
!!***
subroutine op_LaguerreQuadratureRule (nRoots,beta,T,Roots,Weights)

  implicit none

  integer, intent (in) :: nRoots
  real,    intent (in) :: beta
  real,    intent (in) :: T

  real,    intent (out), dimension (1:) :: Roots
  real,    intent (out), dimension (1:) :: Weights

  Roots   = 0.0
  Weights = 0.0

  return
end subroutine op_LaguerreQuadratureRule
