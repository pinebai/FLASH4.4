!!****if* source/physics/materialProperties/Opacity/localAPI/op_LaguerreZeroMoment
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
!! ARGUMENTS
!!
!!  beta    : the exponential for the generalized Laguerre weight function
!!  R       : the upper integration limit (must be > 0)
!!  S       : the exponential factor
!!  MomZero : the value for the 0-th moment
!!
!!***
subroutine op_LaguerreZeroMoment (beta, R, S, MomZero)

  implicit none

  real, intent (in)  :: beta
  real, intent (in)  :: R,S
  real, intent (out) :: MomZero

  MomZero = 0.0

  return
end subroutine op_LaguerreZeroMoment
