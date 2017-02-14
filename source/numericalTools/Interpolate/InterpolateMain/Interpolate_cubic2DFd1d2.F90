!!****if* source/numericalTools/Interpolate/InterpolateMain/Interpolate_cubic2DFd1d2
!!
!! NAME
!!
!!  Interpolate_cubic2DFd1d2
!!
!! SYNOPSIS
!!
!!  Interpolate_cubic2DFd1d2 (real, intent (in) :: a (1:16),
!!                            real, intent (in) :: x,
!!                            real, intent (in) :: y)
!!
!! DESCRIPTION
!!
!!  Calculates the function value and all rescaled 1st and 2nd derivative values for a
!!  pair [x,y] of rescaled [0,1] coordinates and the 16 bicubic expansion coefficients.
!!  The bicubic expansion reads, for one square, in terms of rescaled [0,1] x,y coordinates:
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
!!  The rescaled derivatives are given by the general formula:
!!
!!                 r+s    r  s    3   3                         i-r j-s
!!                d   / dx dy  = sum sum (i) * (j) * a (i,j) * x   y
!!                               i=r j=s    r     s
!!
!!  where the Pochhammer symbols are defined as:
!!
!!                  (i)  = i * (i-1) * (i-2) * ... * (i-r+1)
!!                     r
!!
!!  The rescaled derivatives are therefore also sums of appropriate expansion coefficients
!!  times monomial products. From the general formula we see, that the highest non-zero
!!  derivative is of 6-th order.
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
!!  1) The function is defined as a real array of size 6:
!!
!!           Interpolate_cubic2DFd1d2 (1) = the function value
!!           Interpolate_cubic2DFd1d2 (2) = the rescaled d/dx value
!!           Interpolate_cubic2DFd1d2 (3) = the rescaled d/dy value
!!           Interpolate_cubic2DFd1d2 (4) = the rescaled d2/dx2 value
!!           Interpolate_cubic2DFd1d2 (5) = the rescaled d2/dy2 value
!!           Interpolate_cubic2DFd1d2 (6) = the rescaled d2/dxdy value
!!
!!  2) The code checks, if the supplied pair [x,y] is rescaled.
!!
!!***

function Interpolate_cubic2DFd1d2 (a,x,y)

  use Driver_interface,   ONLY : Driver_abortFlash
  use Logfile_interface,  ONLY : Logfile_stamp

  implicit none

  real, intent (in) :: a (1:16)
  real, intent (in) :: x,y

  logical :: notRescaled

  real    :: x1p5, x2, x6
  real    :: y1p5, y2, y6

  real    :: Interpolate_cubic2DFd1d2 (1:6)  ! declares the function as an array

  real    :: b (1:4)                         ! will hold partial i-index summation values
!
!
!     ...Check status of transmitted rescaled coordinates.
!
!
  notRescaled = (x < 0.0) .or. (y < 0.0) .or. &
                (x > 1.0) .or. (y > 1.0)

  if (notRescaled) then
      call Logfile_stamp     (x, ' = rescaled x coordinate [Interpolate_cubic2DFd1d2]'    )
      call Logfile_stamp     (y, ' = rescaled y coordinate [Interpolate_cubic2DFd1d2]'    )
      call Driver_abortFlash ('[Interpolate_cubic2DFd1d2] ERROR: [x,y] pair not rescaled!')
  end if
!
!
!     ...Calculate needed multiples of the rescaled coordinates.
!
!
  x1p5 = x + (x / 2)
  y1p5 = y + (y / 2)
  x2   = x + x
  y2   = y + y
  x6   = x2 + x2 + x2
  y6   = y2 + y2 + y2
!
!
!     ...Generate the function and derivative values. Intermediate vectors are
!        used to collect common summation contributions. The Horner scheme is
!        applied to calculate the individual polynomial contributions. The
!        multiplication count for each step is listed for information.
!
!                                                                                    ! used for   * count
!
  b (1:4) = ((a (4:16:4) * x + a (3:15:4)) * x + a (2:14:4)) * x + a (1:13:4)        ! f,dy,dy2       12

  Interpolate_cubic2DFd1d2 (1) = ((b (4) * y    + b (3)) * y  + b (2)) * y + b (1)   ! -> f            3
  Interpolate_cubic2DFd1d2 (3) =  (b (4) * y1p5 + b (3)) * y2 + b (2)                ! -> dy           2
  Interpolate_cubic2DFd1d2 (5) =   b (4) * y6   + b (3) + b (3)                      ! -> dy2          1
  
  b (1:4) = (a (4:16:4) * x1p5 + a (3:16:4)) * x2 + a (2:16:4)                       ! dx,dxdy         8

  Interpolate_cubic2DFd1d2 (2) = ((b (4) * y    + b (3)) * y  + b (2)) * y + b (1)   ! -> dx           3
  Interpolate_cubic2DFd1d2 (6) =  (b (4) * y1p5 + b (3)) * y2 + b (2)                ! -> dxdy         2

  b (1:4) = a (4:16:4) * x6 + a (3:16:4) + a (3:16:4)                                ! dx2            16

  Interpolate_cubic2DFd1d2 (4) = ((b (4) * y + b (3)) * y + b (2)) * y + b (1)       ! -> dx2          3
!                                                                                                     --
!                                                                                          total =    50
!
!     ...Ready!
!
!
  return
end function Interpolate_cubic2DFd1d2
