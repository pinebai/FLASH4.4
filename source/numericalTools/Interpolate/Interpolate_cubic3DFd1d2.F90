!!****f* source/numericalTools/Interpolate/Interpolate_cubic3DFd1d2
!!
!! NAME
!!
!!  Interpolate_cubic3DFd1d2
!!
!! SYNOPSIS
!!
!!  Interpolate_cubic3DFd1d2 (real, intent (in) :: a (1:64),
!!                            real, intent (in) :: x,
!!                            real, intent (in) :: y,
!!                            real, intent (in) :: z)
!!
!! DESCRIPTION
!!
!!  Calculates the function value and all rescaled 1st and 2nd derivative values for a
!!  triple [x,y,z] of rescaled [0,1] coordinates and the 64 tricubic expansion coefficients.
!!  The tricubic expansion reads, for one cube, in terms of rescaled [0,1] x,y,z coordinates:
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
!!  The rescaled derivatives are given by the general formula:
!!
!!         r+s+t    r  s  t    3   3   3                                 i-r j-s k-t
!!        d     / dx dy dz  = sum sum sum (i) * (j) * (k) * a (i,j,k) * x   y   z
!!                            i=r j=s k=t    r     s     t
!!
!!  where the Pochhammer symbols are defined as:
!!
!!                  (i)  = i * (i-1) * (i-2) * ... * (i-r+1)
!!                     r
!!
!!  The rescaled derivatives are therefore also sums of appropriate expansion coefficients
!!  times monomial products. From the general formula we see, that the highest non-zero
!!  derivative is of 9-th order.
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
!!  1) The function is defined as a real array of size 10:
!!
!!           Interpolate_cubic3DFd1d2 (1)  = the function value
!!           Interpolate_cubic3DFd1d2 (2)  = the rescaled d/dx value
!!           Interpolate_cubic3DFd1d2 (3)  = the rescaled d/dy value
!!           Interpolate_cubic3DFd1d2 (4)  = the rescaled d/dz value
!!           Interpolate_cubic3DFd1d2 (5)  = the rescaled d2/dx2 value
!!           Interpolate_cubic3DFd1d2 (6)  = the rescaled d2/dy2 value
!!           Interpolate_cubic3DFd1d2 (7)  = the rescaled d2/dz2 value
!!           Interpolate_cubic3DFd1d2 (8)  = the rescaled d2/dxdy value
!!           Interpolate_cubic3DFd1d2 (9)  = the rescaled d2/dxdz value
!!           Interpolate_cubic3DFd1d2 (10) = the rescaled d2/dydz value
!!
!!  2) The code checks, if the supplied triple [x,y,z] is rescaled.
!!
!!***

function Interpolate_cubic3DFd1d2 (a,x,y,z)

  implicit none

  real, intent (in) :: a (1:64)
  real, intent (in) :: x,y,z

  real :: Interpolate_cubic3DFd1d2 (1:10)  ! declares the function as an array

  Interpolate_cubic3DFd1d2 (1:10) = 0.0

  return
end function Interpolate_cubic3DFd1d2
