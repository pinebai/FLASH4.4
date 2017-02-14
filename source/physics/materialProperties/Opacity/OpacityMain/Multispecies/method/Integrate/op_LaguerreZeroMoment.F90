!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Integrate/op_LaguerreZeroMoment
!!
!! NAME
!!
!!  op_LaguerreZeroMoment
!!
!! SYNOPSIS
!!
!!  call op_LaguerreZeroMoment (real (in)  :: beta,
!!                              real (in)  :: R,
!!                              real (in)  :: S,
!!                              real (out) :: MomZero)
!!
!! DESCRIPTION
!!
!!  Condition on the beta value (must be a whole number)
!!  ----------------------------------------------------
!!
!!  The 0-th moment is the integral of x^beta e^(-Sx) over R and 0 and can be evaluated recursively
!!  and exactly for whole number values of beta >= 0. In case beta is positive but not a whole number
!!  one can still apply recursion techniques to bring it down to an integral where beta is bound by
!!  1 > beta > 0, in which case the integral has to be evaluated using approximation methods.
!!  In our case we assume beta is a whole number value and this condition is checked inside the
!!  routine to within the precision accuracy.
!!
!!
!!  Evaluation of the 0-th moment
!!  -----------------------------
!!
!!  The 0-th moment is given by the integral:
!!
!!                                        R
!!                                      /
!!                                     |   beta  -Sx
!!             MomZero (R,S,beta) =    |  x     e    dx
!!                                     |
!!                                    /
!!                                   0
!!
!!  and, for whole number beta's, can be evaluated using the following recursion:
!!
!!                                                                   -RS  b
!!             MomZero (R,S,b) = (1/S) * [b * MomZero (R,S,b-1)  -  e    R ]       b = 1,2,..,beta
!!
!!
!!  starting the recursion with:
!!
!!                                    R
!!                                  /
!!                                 |   -Sx              -RS
!!             MomZero (R,S,0) =   |  e    dx  =  (1 - e   ) / S
!!                                 |
!!                                /
!!                             0
!!
!!  Special care has to be taken for those situations where, due to small exponential exponents,
!!  a significant loss in accuracy can occur when applying the above recursion. This is the case
!!  when the product RS gets small, as will be shown next.
!!
!!  Consider the Taylor series expansion for e^{-Sx} in the above integral. We have:
!!
!!                         e^{-Sx} = 1 - Sx + (Sx)^2/2 - (Sx)^3/6 + ...
!!
!!  Inserting this expansion into the 0-th moment integral, we obtain, after performing the
!!  integration on each term separately the following Taylor-type expansion for the 0-th moment:
!!
!!
!!                                    beta+1    infinity     i    i
!!             MomZero (R,S,beta) =  R            Sum    (-1) (RS)  / [i!(beta+i+1)]
!!                                                i=0
!!
!!
!!  This expression is exact, but only useful for small RS values, such that the sum converges
!!  rapidly. Note, that the absolute values of each term in the summation is actually less than
!!  (RS)^i, because beta >= 0. Hence mantissa accuracy can be easily achieved and the number of
!!  terms that need to be considered can be easily calculated (see the code below).
!!
!!
!! ARGUMENTS
!!
!!  beta    : the exponential for the generalized Laguerre weight function
!!  R       : the upper integration limit (must be > 0)
!!  S       : the exponential factor
!!  MomZero : the value for the 0-th moment
!!
!!***
subroutine op_LaguerreZeroMoment (beta, R, S, MomZero)

  use Driver_interface, ONLY : Driver_abortFlash

  use op_numericsData,  ONLY : zero,one

  implicit none

  real, intent (in)  :: beta
  real, intent (in)  :: R,S
  real, intent (out) :: MomZero

  integer :: b,intbeta
  integer :: i
  integer :: nTerms

  real    :: absRS
  real    :: isNotZeroLimit
  real    :: logR
  real    :: RS
  real    :: Sinv
  real    :: Taylor
  real    :: Term
  real    :: wholeNumberDifference
  real    :: xb,xi

  real, parameter :: smallRSLimit = 0.1
!
!
!   ...Check proper numerical values and size of moment array.
!
!
  if (beta < zero) then
      call Driver_abortFlash ('[op_LaguerreZeroMoment] ERROR: Beta is < 0')
  end if

  if (R <= zero) then
      call Driver_abortFlash ('[op_LaguerreZeroMoment] ERROR: Integration limit R =< 0')
  end if
!
!
!   ...Check whole number consition on beta. Stop, if beta is < 0.
!
!
  intbeta = nint (beta)

  isNotZeroLimit        = epsilon  (one)
  wholeNumberDifference = beta - real (intbeta)

  if (wholeNumberDifference > isNotZeroLimit) then
      call Driver_abortFlash ('[op_LaguerreZeroMoment] ERROR: beta not whole number')
  end if

  if (intbeta < 0) then
      call Driver_abortFlash ('[op_LaguerreZeroMoment] ERROR: beta < 0')
  end if
!
!
!   ...Calculate RS and check, if we can use a Taylor expansion.
!
!
  RS = R * S
  absRS = abs (RS)

  if (absRS < smallRSLimit) then
  
      nTerms = max (1, 1 + int (log (isNotZeroLimit) / log (absRS)))

      Term   = - RS / (beta + one + one)
      Taylor = Term + one / (beta + one)

      do i = 2,nTerms
         xi = real (i)         
         Term   = - RS * Term / xi                      ! accumulate (-1)^i(RS)^i/i! terms
         Taylor = Taylor + Term / (beta + xi + one)
      end do

      MomZero = Taylor * (R ** (beta + 1))

  else
!
!
!   ...Apply the recursion formula.
!
!
      Sinv = one / S

      MomZero = Sinv * (one - exp (-RS))

      if (intbeta > 0) then

          logR = log (R)

          do b = 1,intbeta
             xb = real (b)
             MomZero = Sinv * (xb * MomZero - exp (-RS + xb * logR))
          end do

      end if

  end if
!
!
!   ...Ready! 
!
!
  return
end subroutine op_LaguerreZeroMoment
