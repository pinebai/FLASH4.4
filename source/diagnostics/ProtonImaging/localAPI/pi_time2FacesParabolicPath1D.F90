!!****if* source/diagnostics/ProtonImaging/localAPI/pi_time2FacesParabolicPath1D
!!
!! NAME
!!
!!  pi_time2FacesParabolicPath1D
!!
!! SYNOPSIS
!!
!!  pi_time2FacesParabolicPath1D (real (in) :: pos,
!!                                real (in) :: vel,
!!                                real (in) :: acc,
!!                                real (in) :: minFace,
!!                                real (in) :: maxFace,
!!                                real (in) :: noFaceTime)
!!
!! DESCRIPTION
!!
!!  Given the position, velocity and acceleration of an abstract object and the position
!!  of the minimum and maximum face in 1D, the function calculates the shortest time (excluding
!!  a zero time) the object needs to reach one of the two faces. The underlying path of
!!  the object is parabolic in nature:
!!
!!                     r (t) = r0 + v0 * t + a0 * t^2 / 2
!!
!!  where r0, v0 and a0 are the objects position, velocity and acceleration at time zero and
!!  r (t) is the position of the object at time t. If the object will never hit one of the
!!  two faces, the returned time is the so called 'no face' time, which has to be supplied
!!  by the caller. This enables a time handle for the calling routine in case it needs to
!!  check for successful face hits.
!!
!!  Algorithm:
!!  ---------
!!
!!  In principle it is straightforward to solve the above quadratic equation for time by
!!  inserting the face equations into r (t). The naive approach would be to solve the
!!  two quadratic equations (one for each face) and return the smallest positive root
!!  solutions. This would require 2 sqrt evaluation calls (which are rather expensive)
!!  and an extraction of the smallest positive root out of 4 roots. Since this function
!!  will potentially be called many many times during proton tracing, a more clever way of
!!  doing things results in noticeable performance increase.
!!
!!  The present algorithm analyzes the orientation (positive or negative) of the velocity
!!  and acceleration 'vectors' and decides on which face it will hit first. Thus the
!!  present algorithm needs at most 1 sqrt evaluation and no root checking, at the
!!  expense of a few 'if' decisions.
!!
!! ARGUMENTS
!!
!!  pos        : the object's position
!!  vel        : the object's velocity
!!  acc        : the object's acceleration
!!  minFace    : minimum face position
!!  maxFace    : maximum face position
!!  noFaceTime : the 'no face' time (handle)
!!
!! NOTES
!!        
!!  1) No check is being performed for proper minFace and maxFace values (minFace =< maxFace).
!!
!!  2) The one face case (minFace = maxFace) is allowed and gives the correct hitting time.
!!
!!  3) The position of the object does not have to be between the two faces. It can be
!!     outside of the faces as well.
!!
!!***

real function pi_time2FacesParabolicPath1D (pos, vel, acc, minFace, maxFace, noFaceTime)

  implicit none

  real, intent (in) :: pos, vel, acc
  real, intent (in) :: minFace, maxFace
  real, intent (in) :: noFaceTime

  pi_time2FacesParabolicPath1D = 0.0

  return
end function pi_time2FacesParabolicPath1D
