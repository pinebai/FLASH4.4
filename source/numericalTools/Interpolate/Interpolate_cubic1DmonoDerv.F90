!!****f* source/numericalTools/Interpolate/Interpolate_cubic1DmonoDerv
!!
!! NAME
!!
!!  Interpolate_cubic1DmonoDerv
!!
!! SYNOPSIS
!!
!!  call Interpolate_cubic1DmonoDerv (integer, intent (in)  :: nx,
!!                                    real,    intent (in)  :: f   (0:nx+1),
!!                                    real,    intent (out) :: fx  (1:nx  ))
!!
!! DESCRIPTION
!!
!!  Calculates 1st order x derivatives from a collection of data points, such that
!!  the resulting coefficients for monocubic expansion deliver a monotone 1D line.
!!  While there is always one such monotone line if all derivatives are set equal to
!!  zero, the resulting line is bumpy, especially on data points representing lines of
!!  constant slope. The present routine calculates derivatives, which corrects for this
!!  bumpiness and delivers a more smooth monotonic line. Since the code is dealing with
!!  lines in terms of rescaled [0,1] x coordinates, the resulting derivatives will all
!!  be rescaled.
!!
!!  While there is a whole family of possible monotone lines, the current implementation
!!  is such that it ensures that the cubic expansion for data points alined on a constant
!!  slope will lead to an actual line with that slope. This is achieved by setting the
!!  derivative at a point equal to the minimum of the two differences between the two
!!  neighboring adjacent points and the current point.
!!
!!  Evaluation of the 1st order derivatives needs nearest neighbor data point values
!!  -> we need 1 extra layer of data point values beyond the intended 1D grid.
!!
!! ARGUMENTS
!!
!!  nx        : number of data points in x direction
!!  f   (i)   : data value for i-th grid point
!!  fx  (i)   : 1st order x derivative value for i-th grid point
!!
!! NOTES
!!
!!  ...
!!
!!***

subroutine Interpolate_cubic1DmonoDerv (nx,f,  fx)

  implicit none

  integer, intent (in)  :: nx
  real,    intent (in)  :: f  (0:nx+1)
  real,    intent (out) :: fx (1:nx  )

  fx (1:nx) = 0.0

  return
end subroutine Interpolate_cubic1DmonoDerv
