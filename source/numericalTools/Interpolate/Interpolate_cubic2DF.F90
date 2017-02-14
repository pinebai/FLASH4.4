!!****f* source/numericalTools/Interpolate/Interpolate_cubic2DF
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

  implicit none

  real, intent (in) :: a (1:16)
  real, intent (in) :: x,y

  Interpolate_cubic2DF = 0.0

  return
end function Interpolate_cubic2DF
