!!****if* source/diagnostics/ProtonImaging/localAPI/pi_parabolicPathLength3D
!!
!! NAME
!!
!!  pi_parabolicPathLength3D
!!
!! SYNOPSIS
!!
!!  pi_parabolicPathLength3D (real (in) :: vx,
!!                            real (in) :: vy,
!!                            real (in) :: vz,
!!                            real (in) :: ax,
!!                            real (in) :: ay,
!!                            real (in) :: az,
!!                            real (in) :: T)
!!
!! DESCRIPTION
!!
!!  This function calculates the 3D path length traveled by an object in a certain time
!!  under conditions of constant velocity and acceleration vector. This path length is not
!!  the distance between the initital point and the end point of travel, but rather the
!!  length of the total curve traveled by the object. This path is determined solely by
!!  the velocity and acceleration. The initial position of the object is irrelevant.
!!
!!  The path length L is calculated as an integral:
!!
!!                         T
!!                    L = int |a * t + v| dt
!!                         0
!!
!!  where |a*t+v| is the length of the vector a*t+v at time t. This leads to an integral
!!  involving a square root:
!!
!!                         T
!!                    L = int sqrt (|a|^2 t^2 + 2 a.v t + |v|^2) dt
!!                         0
!!
!!  whose general solution can be stated as:
!!
!!       L = |v|T/2 * [1 + (X + cos)/X * (sqrt[X^2 + 2Xcos + 1] - 1)
!!                       + (1 - cos^2)/X * Ln ((cos + X + sqrt[X^2 + 2Xcos + 1]) / (1 + cos))]
!!
!!  where
!!                   X = |a|T/|v|
!!                 cos = a.v / |a||v|
!!
!!  Several special cases must be considered for computational stability reasons
!!  (E = machine epsilon -> 1.0 + E is different from 1.0):
!!
!!    1) X >= 2/E: In this case the acceleration term |a|T dominates over the
!!                 velocity term |v|, such that the maximum (when cos = 1) and
!!                 minimum (when cos = -1) first order correction of |v|T to
!!                 the distance |a|T^2/2 cannot be accomodated into the mantissa.
!!                 Hence in this case to within computational accuracy:
!!
!!                                L = |a|T^2/2
!!
!!    2) X =< 2E:  Here the velocity term |v| dominates over the acceleration
!!                 term |a|T^2/2, such that the maximum (when cos = 1) and
!!                 minimum (when cos = -1) first order correction of |a|T^2/2 to
!!                 the distance |v|T cannot be accomodated into the mantissa.
!!                 Hence in this case to within computational accuracy:
!!
!!                                L = |v|T
!!
!!    3) cos = 1:  The velocity and acceleration vector are colinear pointing
!!                 in the same direction. The distance is:
!!
!!                                L = |v|T + |a|T^2/2
!!
!!    4) cos = -1: The velocity and acceleration vector are colinear pointing
!!                 in the opposite direction. The distance is:
!!
!!                                L = |v|T - |a|T^2/2
!!
!!  Cases 1) to 4) are checked in that order and if not present, the general formula
!!  is used.
!!
!!  Algorithm:
!!  ---------
!!
!!  In order to avoid costly divisions, the fist two cases are actually checked as:
!!
!!            1) |a|T*E >= 2|v|   and   2) |a|T =< 2E*|v|
!!
!!  Note that evaluation of the distance requires as a minimum the evaluation of
!!  2 square roots (for calculating |v| and |a|) and an additional square root and
!!  natural logarithm evaluation in cases where the general formula must be used.
!!
!! ARGUMENTS
!!
!!  vx : the object's x-component velocity
!!  vy : the object's y-component velocity
!!  vz : the object's z-component velocity
!!  ax : the object's x-component acceleration
!!  ay : the object's y-component acceleration
!!  az : the object's z-component acceleration
!!  T  : the object's travelling time
!!
!! NOTES
!!        
!!  none
!!
!!***

real function pi_parabolicPathLength3D (vx,vy,vz,ax,ay,az,T)

  implicit none

  real, intent (in) :: vx,vy,vz
  real, intent (in) :: ax,ay,az
  real, intent (in) :: T

  pi_parabolicPathLength3D = 0.0

  return
end function pi_parabolicPathLength3D
