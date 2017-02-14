!!****f* source/numericalTools/Interpolate/Interpolate_cubic1DFd1d2
!!
!! NAME
!!
!!  Interpolate_cubic1DFd1d2
!!
!! SYNOPSIS
!!
!!  Interpolate_cubic1DFd1d2 (real, intent (in) :: a (1:4),
!!                            real, intent (in) :: x)
!!
!! DESCRIPTION
!!
!!  Calculates the function value and the rescaled 1st and 2nd derivative value for a
!!  single [x] rescaled [0,1] coordinate and the 4 monocubic expansion coefficients.
!!  The monocubic expansion reads, for one line, in terms of the rescaled [0,1] x coordinate:
!!
!!                                 3          i
!!                        F (x) = sum  a (i) x
!!                                i=0
!!
!!  The location index of the a (i) inside the 4-dimensional vector is:
!!
!!                location index of a (i)  =  1 + i
!!
!!  The rescaled derivatives are given by the general formula:
!!
!!                 r    r    3                  i-r
!!                d / dx  = sum (i)  * a (i) * x
!!                          i=r    r
!!
!!  where the Pochhammer symbols are defined as:
!!
!!                  (i)  = i * (i-1) * (i-2) * ... * (i-r+1)
!!                     r
!!
!!  From the derivative formula we see, that the highest non-zero derivative is of
!!  3-rd order.
!!
!! ARGUMENTS
!!
!!  a (i) : the i-th monocubic expansion coefficient
!!  x     : rescaled [0,1] x coordinate
!!
!! NOTES
!!
!!  1) The function is defined as a real array of size 3:
!!
!!           Interpolate_cubic1DFd1d2 (1) = the function value
!!           Interpolate_cubic1DFd1d2 (2) = the rescaled d/dx value
!!           Interpolate_cubic1DFd1d2 (3) = the rescaled d2/dx2 value
!!
!!  2) The code checks, if the supplied coordinate [x] is rescaled.
!!
!!***

function Interpolate_cubic1DFd1d2 (a,x)

  implicit none

  real, intent (in) :: a (1:4)
  real, intent (in) :: x

  real :: Interpolate_cubic1DFd1d2 (1:3)  ! declares the function as an array

  Interpolate_cubic1DFd1d2 (1:3) = 0.0

  return
end function Interpolate_cubic1DFd1d2
