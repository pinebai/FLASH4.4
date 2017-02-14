!!****if* source/numericalTools/Interpolate/InterpolateMain/Interpolate_cubic2DF
!!
!! NAME
!!
!!  Interpolate_cubic2DF
!!
!! SYNOPSIS
!!
!!  Interpolate_cubic2DF (real, intent (in) :: a (1:16),
!!                        real, intent (in) :: x,
!!                        real, intent (in) :: y)
!!
!! DESCRIPTION
!!
!!  Calculates the function value for a pair [x,y] of rescaled [0,1] coordinates and the
!!  16 bicubic expansion coefficients. The bicubic expansion reads, for one square, in
!!  terms of rescaled [0,1] x,y coordinates:
!!
!!                                   3   3            i j
!!                        F (x,y) = sum sum  a (i,j) x y
!!                                  i=0 j=0
!!
!!  The order of the supplied expansion coefficients a (i,j) must be such, that the
!!  j index has the highest ranking, followed by the i index. The overall location index
!!  of the a (i,j) inside the 16-dimensional vector is given by the following formula:
!!
!!                location index of (i,j)  =  1 + i + 4j
!!
!!  Since this function is (potentially) called many times from external applications,
!!  efficiency is key here and intermediate common summation terms are reused as much as
!!  possible. The strategy is partial summation and reduction at each index summation stage.
!!  The individual x- and y-coordinate cubic polynomial sections are always evaluated
!!  using the Horner scheme to minimize accumulation of computation rounding errors.
!!
!! ARGUMENTS
!!
!!  a (i) : the i-th bicubic expansion coefficient
!!  x     : rescaled [0,1] x coordinate
!!  y     : rescaled [0,1] y coordinate
!!
!! NOTES
!!
!!  1) The code checks, if the supplied pair [x,y] is rescaled.
!!
!!***

real function Interpolate_cubic2DF (a,x,y)

  use Driver_interface,   ONLY : Driver_abortFlash
  use Logfile_interface,  ONLY : Logfile_stamp

  implicit none

  real, intent (in) :: a (1:16)
  real, intent (in) :: x,y

  logical :: notRescaled

  real    :: b (1:4)       ! will hold partial i-index  summation values
!
!
!     ...Check status of transmitted rescaled coordinates.
!
!
  notRescaled = (x < 0.0) .or. (y < 0.0) .or. &
                (x > 1.0) .or. (y > 1.0)

  if (notRescaled) then
      call Logfile_stamp     (x, ' = rescaled x coordinate [Interpolate_cubic2DF]'    )
      call Logfile_stamp     (y, ' = rescaled y coordinate [Interpolate_cubic2DF]'    )
      call Driver_abortFlash ('[Interpolate_cubic2DF] ERROR: [x,y] pair not rescaled!')
  end if
!
!
!     ...Generate the function value. Intermediate vectors are used to collect common
!        summation contributions. The Horner scheme is applied to calculate the individual
!        polynomial contributions. The multiplication count for each step is listed for
!        information.
!
!                                                                               ! used for   * count
!
  b (1:4) = ((a (4:16:4) * x + a (3:15:4)) * x + a (2:14:4)) * x + a (1:13:4)   ! f              12

  Interpolate_cubic2DF = ((b (4) * y + b (3)) * y + b (2)) * y + b (1)          ! -> f            3
!                                                                                                --
!                                                                                     total =    15
!
!     ...Ready!
!
!
  return
end function Interpolate_cubic2DF
