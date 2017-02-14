!!****f* source/numericalTools/Interpolate/Interpolate_cubic3DF
!!
!! NAME
!!
!!  Interpolate_cubic3DF
!!
!! SYNOPSIS
!!
!!  Interpolate_cubic3DF (real, intent (in) :: a (1:64),
!!                        real, intent (in) :: x,
!!                        real, intent (in) :: y,
!!                        real, intent (in) :: z)
!!
!! DESCRIPTION
!!
!!  Calculates the function value for a triple [x,y,z] of rescaled [0,1] coordinates and the
!!  64 tricubic expansion coefficients. The tricubic expansion reads, for one cube, in terms
!!  of rescaled [0,1] x,y,z coordinates:
!!
!!                                   3   3   3              i j k
!!                      F (x,y,z) = sum sum sum  a (i,j,k) x y z
!!                                  i=0 j=0 k=0
!!
!!  The order of the supplied expansion coefficients a (i,j,k) must be such, that the
!!  k index has the highest ranking, followed by the j index and the i index. The overall
!!  location index of the a (i,j,k) inside the 64-dimensional vector is given by the
!!  following formula:
!!
!!                location index of (i,j,k)  =  1 + i + 4j + 16k
!!
!!  Since this function is (potentially) called many times from external applications,
!!  efficiency is key here and intermediate common summation terms are reused as much as
!!  possible. The strategy is partial summation and reduction at each index summation stage.
!!  The individual x-,y- and z-coordinate cubic polynomial sections are always evaluated
!!  using the Horner scheme to minimize accumulation of computation rounding errors.
!!
!! ARGUMENTS
!!
!!  a (i) : the i-th tricubic expansion coefficient
!!  x     : rescaled [0,1] x coordinate
!!  y     : rescaled [0,1] y coordinate
!!  z     : rescaled [0,1] z coordinate
!!
!! NOTES
!!
!!  The code checks, if the supplied triple [x,y,z] is rescaled.
!!
!!***

real function Interpolate_cubic3DF (a,x,y,z)

  implicit none

  real, intent (in) :: a (1:64)
  real, intent (in) :: x,y,z

  Interpolate_cubic3DF = 0.0

  return
end function Interpolate_cubic3DF
