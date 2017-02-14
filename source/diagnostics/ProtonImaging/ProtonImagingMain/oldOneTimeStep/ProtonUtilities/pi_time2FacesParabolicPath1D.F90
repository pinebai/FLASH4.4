!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/oldOneTimeStep/ProtonUtilities/pi_time2FacesParabolicPath1D
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

  real :: b,c,t
!
!
!     ...Determine the time to reach one of the two faces.
!
!
  t = noFaceTime                                   ! set the default time

  if (vel > 0.0) then

      if (acc > 0.0) then                          ! velocity and acceleration in same direction

          if (minFace > pos) then                  ! no equal sign (exclude zero times)
 
              c = 2 * (minFace - pos)              ! + ve
              b = vel * vel + acc * c              ! + ve
              t = c / (vel + sqrt (b))             ! the smaller positive time (1st minFace hit)

          else if (maxFace > pos) then             ! no equal sign (exclude zero times)

              c = 2 * (maxFace - pos)              ! + ve
              b = vel * vel + acc * c              ! + ve
              t = c / (vel + sqrt (b))             ! the smaller positive time (1st maxFace hit)

          end if

      else if (acc < 0.0) then                     ! velocity and acceleration in opposite direction

          if (minFace > pos) then                  ! no equal sign (exclude zero times)

              c = 2 * (minFace - pos)              ! + ve
              b = vel * vel + acc * c              ! +/- ve
              if (b >= 0.0) then                   ! it hits minFace
                  t = c / (vel + sqrt (b))         ! the smaller positive time (1st minFace hit)
              end if

          else if (maxFace > pos) then             ! no equal sign (exclude zero times)

              c = 2 * (maxFace - pos)              ! + ve
              b = vel * vel + acc * c              ! +/- ve
              if (b >= 0.0) then                   ! it hits maxFace
                  t = c / (vel + sqrt (b))         ! the smaller positive time (1st maxFace hit)
              else                                 ! it must hit minFace
                  c = 2 * (minFace - pos)          ! - ve
                  b = vel * vel + acc * c          ! + ve
                  t = - (vel + sqrt (b)) / acc     ! the only positive time (1st minFace hit)
              end if

          else if (maxFace < pos) then             ! no equal sign (exclude zero times)

              c = 2 * (maxFace - pos)              ! - ve
              b = vel * vel + acc * c              ! + ve
              t = - (vel + sqrt (b)) / acc         ! the only positive time (1st maxFace hit)

          end if

      else                                         ! no acceleration (solve linear time equation)

          if (minFace > pos) then                  ! no equal sign (exclude zero times)
              t = (minFace - pos) / vel            ! the only positive time (1st minFace hit)
          else if (maxFace > pos) then             ! no equal sign (exclude zero times)
              t = (maxFace - pos) / vel            ! the only positive time (1st maxFace hit)
          end if

      end if

  else if (vel < 0.0) then

      if (acc < 0.0) then                          ! velocity and acceleration in same direction

          if (maxFace < pos) then                  ! no equal sign (exclude zero times)

              c = 2 * (maxFace - pos)              ! - ve
              b = vel * vel + acc * c              ! + ve
              t = c / (vel - sqrt (b))             ! the smaller positive time (1st maxFace hit)

          else if (minFace < pos) then             ! no equal sign (exclude zero times)

              c = 2 * (minFace - pos)              ! - ve
              b = vel * vel + acc * c              ! + ve
              t = c / (vel - sqrt (b))             ! the smaller positive time (1st minFace hit)

          end if

      else if (acc > 0.0) then                     ! velocity and acceleration in opposite direction

          if (maxFace < pos) then                  ! no equal sign (exclude zero times)

              c = 2 * (maxFace - pos)              ! - ve
              b = vel * vel + acc * c              ! +/- ve
              if (b >= 0.0) then                   ! it hits maxFace
                  t = c / (vel - sqrt (b))         ! the smaller positive time (1st maxFace hit)
              end if

          else if (minFace < pos) then             ! no equal sign (exclude zero times)

              c = 2 * (minFace - pos)              ! - ve
              b = vel * vel + acc * c              ! +/- ve
              if (b >= 0.0) then                   ! it hits minFace
                  t = c / (vel - sqrt (b))         ! the smaller positive time (1st minFace hit)
              else                                 ! it must hit maxFace
                  c = 2 * (maxFace - pos)          ! + ve
                  b = vel * vel + acc * c          ! + ve
                  t = (- vel + sqrt (b)) / acc     ! the only positive time (1st maxFace hit)
              end if

          else if (minFace > pos) then             ! no equal sign (exclude zero times)

              c = 2 * (minFace - pos)              ! + ve
              b = vel * vel + acc * c              ! + ve
              t = (- vel + sqrt (b)) / acc         ! the only positive time (1st minFace hit)

          end if

      else                                         ! no acceleration (solve linear time equation)

          if (maxFace < pos) then                  ! no equal sign (exclude zero times)
              t = (maxFace - pos) / vel            ! the only positive time (1st maxFace hit)
          else if (minFace < pos) then             ! no equal sign (exclude zero times)
              t = (minFace - pos) / vel            ! the only positive time (1st minFace hit)
          end if

      end if

  else

      if (acc > 0.0) then                          ! only positive acceleration (no velocity)

          if (minFace > pos) then                  ! no equal sign (exclude zero times)

              t = (minFace - pos) / acc            ! + ve
              t = sqrt (t + t)                     ! the time to hit minFace

          else if (maxFace > pos) then             ! no equal sign (exclude zero times)

              t = (maxFace - pos) / acc            ! + ve
              t = sqrt (t + t)                     ! the time to hit maxFace

          end if

      else if (acc < 0.0) then                     ! only negative acceleration (no velocity)
      
          if (maxFace < pos) then                  ! no equal sign (exclude zero times)

              t = (maxFace - pos) / acc            ! + ve
              t = sqrt (t + t)                     ! the time to hit maxFace

          else if (minFace < pos) then             ! no equal sign (exclude zero times)

              t = (minFace - pos) / acc            ! + ve
              t = sqrt (t + t)                     ! the time to hit minFace

          end if

      end if

  end if

  pi_time2FacesParabolicPath1D = t
!
!
!     ...Ready!
!
!
  return
end function pi_time2FacesParabolicPath1D
