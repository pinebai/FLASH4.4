!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Integrate/op_LaguerreQuadratureRule
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

  use Driver_interface, ONLY : Driver_abortFlash

  use op_interface,     ONLY : op_LaguerreMoments,       &
                               op_LaguerreZeroMoment,    &
                               op_LaguerreCoefficients,  &
                               op_shJacobi3TermMoments,  &
                               op_shJacobiCoefficients,  &
                               op_OrthoPolyCoefficients, &
                               op_generateRootsWeights,  &
                               op_writeQuadratureData

  use op_integrateData, ONLY : op_maxRoots,              &
                               op_maxWork,               &
                               op_Moments,               &
                               op_JmatrixDiagonals,      &
                               op_JmatrixOffdiagonals,   &
                               op_AuxPolynomialA,        &
                               op_AuxPolynomialB,        &
                               op_OrthoPolynomialA,      &
                               op_OrthoPolynomialB,      &
                               op_work1,                 &
                               op_work2,                 &
                               op_work3,                 &
                               op_printQuadratureData

  use op_numericsData,  ONLY : zero,one

  implicit none

  integer, intent (in)  :: nRoots
  real,    intent (in)  :: beta
  real,    intent (in)  :: T

  real,    intent (out) :: Roots   (1:nRoots)
  real,    intent (out) :: Weights (1:nRoots)

  logical :: Ttiny
  logical :: Tsmall
  logical :: Tlarge
  logical :: Thuge

  integer :: nMoments
  integer :: orthogonalPolynomial

  integer, parameter :: known    = 1
  integer, parameter :: unknown  = 0

  real :: LimitTsmall
  real :: MomZero
  real :: p,q
  real :: Tfraction

  real, parameter :: LimitTtiny   =   1.E-12
  real, parameter :: LimitTlarge  = 150.0
!
!
!   ...Check, if the # of roots wanted conformes with the array space allocated.
!      Currently, only up to 20 roots can be evaluated safely with high precision.
!
!
  if (nRoots > op_maxRoots) then
      call Driver_abortFlash ('[op_LaguerreQuadratureRule] ERROR: # of roots too large')
  end if

  if (nRoots > 40) then
      call Driver_abortFlash ('[op_LaguerreQuadratureRule] ERROR: # of roots limited to =< 40')
  end if
!
!
!   ...Check numerical values passed in argument.
!
!
  if (nRoots < 1) then
      call Driver_abortFlash ('[op_LaguerreQuadratureRule] ERROR: invalid # of roots')
  end if

  if (beta < 0) then
      call Driver_abortFlash ('[op_LaguerreQuadratureRule] ERROR: Beta is < 0')
  end if

  if (T <= zero) then
      call Driver_abortFlash ('[op_LaguerreQuadratureRule] ERROR: Integration limit T =< 0')
  end if

  if (size (Roots) < nRoots) then
      call Driver_abortFlash ('[op_LaguerreQuadratureRule] ERROR: Roots array size too small')
  end if

  if (size (Weights) < nRoots) then
      call Driver_abortFlash ('[op_LaguerreQuadratureRule] ERROR: Weights array size too small')
  end if
!
!
!   ...Initialize some variables (even if they are not needed!).
!
!
  p = zero
  q = zero

  Ttiny  = .false.
  Tsmall = .false.
  Tlarge = .false.
  Thuge  = .false.

  if (nRoots <= 20) then
      LimitTsmall  =  60.0
  else if (nRoots <= 30) then
      LimitTsmall  =  100.0
  else if (nRoots <= 40) then
      LimitTsmall  =  150.0
  end if

  nMoments = 2 * nRoots - 1
!
!
!   ...Calculate the modified moments and the 3-term recursion coefficients of the
!      auxilliary classical polynomials. Different appraoches are used depending on
!      the size of the upper integration limit T.
!
!
  if (T < LimitTtiny) then

      Ttiny = .true.
      p = beta + one

      call op_LaguerreZeroMoment    (beta,          &
                                     one,T,         &
                                             MomZero)

      call op_shJacobiCoefficients  (nRoots,                      &
                                     p,p,                         &
                                             op_OrthoPolynomialA, &
                                             op_OrthoPolynomialB  )

      orthogonalPolynomial = known

  else if (T < LimitTsmall) then

      Tsmall = .true.
      Tfraction = one - (one / LimitTsmall) * real (nRoots)    ! this holds for nRoots =< 40
 
      p = (T * Tfraction) + beta + one
      q = beta + one

      write (*,*) ' p = ',p

      call op_shJacobi3TermMoments  (nMoments,              &
                                     op_maxWork,            &
                                     p,                     &
                                     beta,                  &
                                     T,                     &
                                     op_work1,              &
                                     op_work2,              &
                                     op_work3,              &
                                                 MomZero,   &
                                                 op_Moments )

      call op_shJacobiCoefficients  (nMoments,                      &
                                     p,q,                           &
                                                 op_AuxPolynomialA, &
                                                 op_AuxPolynomialB  )

      orthogonalPolynomial = unknown

  else if (T < LimitTlarge) then

      Tlarge = .true.

      call op_LaguerreMoments       (nMoments,              &
                                     op_maxWork,            &
                                     beta,                  &
                                     beta,                  &
                                     T,                     &
                                     op_work1,              &
                                                 MomZero,   &
                                                 op_Moments )

      call op_LaguerreCoefficients  (nMoments,                      &
                                     beta,                          &
                                                 op_AuxPolynomialA, &
                                                 op_AuxPolynomialB  )

      orthogonalPolynomial = unknown

  else

      Thuge = .true.

      call op_LaguerreZeroMoment    (beta,          &
                                     T,one,         &
                                             MomZero)

      call op_LaguerreCoefficients  (nRoots,                      &
                                     beta,                        &
                                             op_OrthoPolynomialA, &
                                             op_OrthoPolynomialB  )

      orthogonalPolynomial = known

  end if
!
!
!   ...Generate the 3-term recursion coefficients for the orthogonal polynomials
!      (if not known at this stage).
!
!
  if (orthogonalPolynomial == unknown) then

      call op_OrthoPolyCoefficients (nRoots,                                &
                                     nMoments,                              &
                                     op_maxWork,                            &
                                     op_Moments,                            &
                                     op_AuxPolynomialA,                     &
                                     op_AuxPolynomialB,                     &
                                     op_work1,                              &
                                     op_work2,                              &
                                                       op_OrthoPolynomialA, &
                                                       op_OrthoPolynomialB  )
  end if
!
!
!   ...Set up the terminal Jacobi matrix from the set of obtained orthogonal
!      polynomial recursion coefficients and generate the roots and weights.
!      If the shifted Jacobi polynomials were used, do not forget to rescale
!      the roots and weights.
!
!
  op_JmatrixDiagonals    (1:nRoots  ) =       op_OrthoPolynomialA (1:nRoots)
  op_JmatrixOffdiagonals (1:nRoots-1) = sqrt (op_OrthoPolynomialB (2:nRoots))

  call op_generateRootsWeights (nRoots,                               &
                                op_maxWork,                           &
                                MomZero,                              &
                                op_work1,                             &
                                op_JmatrixDiagonals,                  &
                                op_JmatrixOffdiagonals,               &
                                                        Roots,Weights )

  if (Ttiny .or. Tsmall) then
      Roots   (1:nRoots) =                 T * Roots   (1:nRoots)
      Weights (1:nRoots) = T ** (beta + one) * Weights (1:nRoots)
  end if
!
!
!   ...If requested, print out the quadrature data.
!
!
  if (op_printQuadratureData) then

      call op_writeQuadratureData (nRoots,   &
                                   nMoments, &
                                   Ttiny,    &
                                   Tsmall,   &
                                   Tlarge,   &
                                   Thuge,    &
                                   T,        &
                                   p,q,      &
                                   beta,     &
                                   MomZero,  &
                                   Roots,    &
                                   Weights   )
  end if
!
!
!   ...Ready! 
!
!
  return
end subroutine op_LaguerreQuadratureRule
