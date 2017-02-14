!!****f* source/numericalTools/Interpolate/Interpolate_cubic2DmonoDerv
!!
!! NAME
!!
!!  Interpolate_cubic2DmonoDerv
!!
!! SYNOPSIS
!!
!!  call Interpolate_cubic2DmonoDerv (integer, intent (in)  :: nx,
!!                                    integer, intent (in)  :: ny,
!!                                    real,    intent (in)  :: f   (-1:nx+2,-1:ny+2),
!!                                    real,    intent (out) :: fx  ( 0:nx+1, 0:ny+1),
!!                                    real,    intent (out) :: fy  ( 0:nx+1, 0:ny+1),
!!                                    real,    intent (out) :: fxy ( 1:nx  , 1:ny  ))
!!
!! DESCRIPTION
!!
!!  Calculates 1st order x and y derivatives and 2nd order xy mixed derivatives from
!!  a collection of square data points, such that the resulting coefficients for bicubic
!!  expansion deliver a monotone 2D surface. While there is always one such monotone
!!  surface if all derivatives are set equal to zero, the resulting 2D surface is bumpy,
!!  especially on data points representing surfaces of constant slope in either x or
!!  y direction. The present routine calculates derivatives, which corrects for this
!!  bumpiness and delivers a more smooth monotonic surface. Since the code is dealing
!!  with squares in terms of rescaled [0,1] x,y coordinates, the resulting derivatives will
!!  all be rescaled.
!!
!!  Labeling convention of all variables according to square grid point location:
!!
!!
!!                                              y (j index)
!!
!!                                              |
!!                                              |
!!
!!                                 mp --------- 0p -------- pp
!!
!!                                 |            |            |
!!                                 |            |            |
!!                                 |            |            |
!!                                 |            |            |
!!
!!                                 m0 --------- 00 -------- p0 -----> x (i index)
!!
!!                                 |            |            |
!!                                 |            |            |
!!                                 |            |            |
!!                                 |            |            |
!!
!!                                 mm --------- 0m -------- pm
!!
!!
!!  For evaluating the 2nd order mixed derivative at a point (i,j), the nearest neighbor
!!  1st order x and y derivatives are needed -> we need an extra layer of 1st order derivatives
!!  beyond the intended 2D grid. Evaluation of the 1st order derivatives needs next neighbor
!!  data point values -> we need 2 extra layers of data point values beyond the intended 2D grid.
!!
!! ARGUMENTS
!!
!!  nx        : number of data points in x direction
!!  ny        : number of data points in y direction
!!  f   (i,j) : data value for i,j-th grid point
!!  fx  (i,j) : 1st order x derivative value for i,j-th grid point
!!  fy  (i,j) : 1st order y derivative value for i,j-th grid point
!!  fxy (i,j) : 2nd order mixed xy derivative value for i,j-th grid point
!!
!! NOTES
!!
!!  ...
!!
!!***

subroutine Interpolate_cubic2DmonoDerv (nx, ny,        &
                                        f,             &
                                                   fx, &
                                                   fy, &
                                                   fxy )

  implicit none

  integer, intent (in)  :: nx, ny
  real,    intent (in)  :: f   (-1:nx+2 , -1:ny+2)
  real,    intent (out) :: fx  ( 0:nx+1 ,  0:ny+1)
  real,    intent (out) :: fy  ( 0:nx+1 ,  0:ny+1)
  real,    intent (out) :: fxy ( 1:nx   ,  1:ny  )

  fx  (0:nx+1 , 0:ny+1) = 0.0
  fy  (0:nx+1 , 0:ny+1) = 0.0
  fxy (1:nx   , 1:ny  ) = 0.0

  return
end subroutine Interpolate_cubic2DmonoDerv
