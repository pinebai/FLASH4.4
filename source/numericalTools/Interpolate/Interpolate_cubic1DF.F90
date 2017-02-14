!!****f* source/numericalTools/Interpolate/Interpolate_cubic1DF
!!
!! NAME
!!
!!  Interpolate_cubic1DF
!!
!! SYNOPSIS
!!
!!  Interpolate_cubic1DF (real, intent (in) :: a (1:4),
!!                        real, intent (in) :: x)
!!
!! DESCRIPTION
!!
!!  Calculates the function value for a single [x] rescaled [0,1] coordinate and the
!!  4 monocubic expansion coefficients. The monocubic expansion reads, for one line,
!!  in terms of the rescaled [0,1] x coordinate:
!!
!!                                 3          i
!!                        F (x) = sum  a (i) x
!!                                i=0
!!
!!  The location index of the a (i) inside the 4-dimensional vector is:
!!
!!                location index of a (i)  =  1 + i
!!
!! ARGUMENTS
!!
!!  a (i) : the i-th monocubic expansion coefficient
!!  x     : rescaled [0,1] x coordinate
!!
!! NOTES
!!
!!  1) The code checks, if the supplied coordinate [x] is rescaled.
!!
!!***

real function Interpolate_cubic1DF (a,x)

  implicit none

  real, intent (in) :: a (1:4)
  real, intent (in) :: x

  Interpolate_cubic1DF = 0.0

  return
end function Interpolate_cubic1DF
