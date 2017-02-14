!!****f* source/numericalTools/Interpolate/Interpolate_cubic2DFd1d2
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

  implicit none

  real, intent (in) :: a (1:16)
  real, intent (in) :: x,y

  real :: Interpolate_cubic2DFd1d2 (1:6)  ! declares the function as an array

  Interpolate_cubic2DFd1d2 (1:6) = 0.0

  return
end function Interpolate_cubic2DFd1d2
