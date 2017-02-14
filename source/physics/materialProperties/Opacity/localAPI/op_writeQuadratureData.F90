!!****if* source/physics/materialProperties/Opacity/localAPI/op_writeQuadratureData
!!
!! NAME
!!
!!  op_writeQuadratureData
!!
!! SYNOPSIS
!!
!!  call op_writeQuadratureData (integer (in) :: nRoots,
!!                               integer (in) :: nMoments,
!!                               logical (in) :: Ttiny,
!!                               logical (in) :: Tsmall,
!!                               logical (in) :: Tlarge,
!!                               logical (in) :: Thuge,
!!                               real    (in) :: T,
!!                               real    (in) :: p,
!!                               real    (in) :: q,
!!                               real    (in) :: beta,
!!                               real    (in) :: MomZero,
!!                               real    (in) :: Roots   (1:nRoots),
!!                               real    (in) :: Weights (1:nRoots))
!!
!! DESCRIPTION
!!
!!  Utility routine, which prints detailed info regarding the presently calculated
!!  quadrature rule to a text file. The information is written out to a file
!!  named basenm_opacityQuadratureData.txt, where basenm is the runtime parameter for
!!  output file names. The file is appended at each time for each quadrature.
!!  The routine is mainly meant for checking purpose only.
!!
!! ARGUMENTS
!!
!!  nRoots   : the number of roots and weights calculated
!!  nMoments : the number of moments calculated
!!  Ttiny    : indicator, if the integration limit T was considered tiny
!!  Tsmall   : indicator, if the integration limit T was considered small
!!  Tlarge   : indicator, if the integration limit T was considered large
!!  Thuge    : indicator, if the integration limit T was considered huge
!!  T        : the upper integration limit
!!  p        : the p parameter for the shifted Jacobi polynomials G(p,q,x)
!!  q        : the q parameter for the shifted Jacobi polynomials G(p,q,x)
!!  beta     : the exponent in x^beta for the Laguerre-type weight
!!  MomZero  : the 0-th moment
!!  Roots    : the roots for i = 1,2,...,nRoots
!!  Weights  : the weights for i = 1,2,...,nRoots
!!
!!***
subroutine op_writeQuadratureData (nRoots, nMoments,                  &
                                   Ttiny, Tsmall, Tlarge, Thuge, T,   &
                                   p, q, beta, MomZero, Roots, Weights)

  implicit none

  integer, intent (in) :: nRoots
  integer, intent (in) :: nMoments
  logical, intent (in) :: Ttiny
  logical, intent (in) :: Tsmall
  logical, intent (in) :: Tlarge
  logical, intent (in) :: Thuge
  real,    intent (in) :: T
  real,    intent (in) :: p,q
  real,    intent (in) :: beta
  real,    intent (in) :: MomZero

  real,    intent (in) :: Roots   (1:nRoots)
  real,    intent (in) :: Weights (1:nRoots)

  return
end subroutine op_writeQuadratureData
