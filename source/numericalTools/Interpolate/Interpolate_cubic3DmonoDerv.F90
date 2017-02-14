!!****f* source/numericalTools/Interpolate/Interpolate_cubic3DmonoDerv
!!
!! NAME
!!
!!  Interpolate_cubic3DmonoDerv
!!
!! SYNOPSIS
!!
!!  call Interpolate_cubic3DmonoDerv (integer, intent (in)  :: nx,
!!                                    integer, intent (in)  :: ny,
!!                                    integer, intent (in)  :: nz,
!!                                    real,    intent (in)  :: f    (-2:nx+3 , -2:ny+3 , -2:nz+3),
!!                                    real,    intent (out) :: fx   (-1:nx+2 , -1:ny+2 , -1:nz+2),
!!                                    real,    intent (out) :: fy   (-1:nx+2 , -1:ny+2 , -1:nz+2),
!!                                    real,    intent (out) :: fz   (-1:nx+2 , -1:ny+2 , -1:nz+2),
!!                                    real,    intent (out) :: fxy  ( 0:nx+1 ,  0:ny+1 ,  0:nz+1),
!!                                    real,    intent (out) :: fxz  ( 0:nx+1 ,  0:ny+1 ,  0:nz+1),
!!                                    real,    intent (out) :: fyz  ( 0:nx+1 ,  0:ny+1 ,  0:nz+1),
!!                                    real,    intent (out) :: fxyz ( 1:nx   ,  1:ny   ,  1:nz  ))
!!
!! DESCRIPTION
!!
!!  Calculates 1st order x,y,z derivatives, 2nd order xy,xz,yz mixed derivatives and
!!  3rd order xyz mixed derivatives from a collection of cube data points, such that the
!!  resulting coefficients for a triicubic expansion deliver a monotone 3D surface. While
!!  there is always one such monotone surface if all derivatives are set equal to zero, the
!!  resulting 3D surface is bumpy, especially on data points representing surfaces of constant
!!  slope in either x,y or z direction. The present routine calculates derivatives, which
!!  corrects for this bumpiness and delivers a more smooth monotonic surface. Since the code
!!  is dealing with cubes in terms of rescaled [0,1] x,y,z coordinates, the resulting derivatives
!!  will all be rescaled.
!!
!!  Labeling convention of all variables according to cube grid point location (only the
!!  positive y-axis region is shown for clarity):
!!
!!
!!
!!                           mpp -------- 0pp -------- ppp
!!                           /|           /|           /|
!!                          / |          / |          / |
!!                         /  |         /  |         /  |
!!                        /   |        /   |        /   |
!!                       /    |       /    |       /    |        z (k index)
!!                     m0p -------- 00p -------- p0p    |
!!                      |    mp0 ----|--- 0p0 ----|--- pp0           |
!!                      |    /|      |    /|      |    /|            |  y (j index)
!!                      |   / |      |   / |      |   / |            |
!!                      |  /  |      |  /  |      |  /  |            |  /
!!                      | /   |      | /   |      | /   |            | /
!!                      |/    |      |/    |      |/    |            |/
!!                     m00 -------- 000 -------- p00    |             --------- x (i index)
!!                      |    mpm ----|--- 0pm ----|--- ppm
!!                      |    /       |    /       |    /
!!                      |   /        |   /        |   /
!!                      |  /         |  /         |  /
!!                      | /          | /          | /
!!                      |/           |/           |/
!!                     m0m -------- 00m -------- p0m
!!
!!
!!
!!
!!  For evaluating the 3rd order mixed derivative at a point (i,j,k), the nearest neighbor
!!  2nd order derivatives are needed -> we need an extra layer of 2nd order derivatives
!!  beyond the intended 3D grid. Evaluation of the 2nd order derivatives needs nearest neighbor
!!  1st order derivatives -> we need 2 extra layers of 1st order derivatives beyond the
!!  intended 3D grid. Evaluation of the 1st order derivatives needs nearest neighbor
!!  data point values -> we need 3 extra layers of data point values beyond the intended
!!  3D grid.
!!
!! ARGUMENTS
!!
!!  nx           : number of intended 3D grid points in x direction
!!  ny           : number of intended 3D grid points in y direction
!!  nz           : number of intended 3D grid points in z direction
!!  f    (i,j,k) : data value for i,j,k-th grid point
!!  fx   (i,j,k) : 1st order x derivative value for i,j,k-th grid point
!!  fy   (i,j,k) : 1st order y derivative value for i,j,k-th grid point
!!  fz   (i,j,k) : 1st order z derivative value for i,j,k-th grid point
!!  fxy  (i,j,k) : 2nd order mixed xy derivative value for i,j,k-th grid point
!!  fxz  (i,j,k) : 2nd order mixed xz derivative value for i,j,k-th grid point
!!  fyz  (i,j,k) : 2nd order mixed yz derivative value for i,j,k-th grid point
!!  fxyz (i,j,k) : 3rd order mixed xyz derivative value for i,j,k-th grid point
!!
!! NOTES
!!
!!  ...
!!
!!***

subroutine Interpolate_cubic3DmonoDerv (nx, ny, nz,           &
                                        f,                    &
                                                 fx,fy,fz,    &
                                                 fxy,fxz,fyz, &
                                                 fxyz         )

  implicit none

  integer, intent (in)  :: nx, ny, nz
  real,    intent (in)  :: f    (-2:nx+3 , -2:ny+3 , -2:nz+3)
  real,    intent (out) :: fx   (-1:nx+2 , -1:ny+2 , -1:nz+2)
  real,    intent (out) :: fy   (-1:nx+2 , -1:ny+2 , -1:nz+2)
  real,    intent (out) :: fz   (-1:nx+2 , -1:ny+2 , -1:nz+2)
  real,    intent (out) :: fxy  ( 0:nx+1 ,  0:ny+1 ,  0:nz+1)
  real,    intent (out) :: fxz  ( 0:nx+1 ,  0:ny+1 ,  0:nz+1)
  real,    intent (out) :: fyz  ( 0:nx+1 ,  0:ny+1 ,  0:nz+1)
  real,    intent (out) :: fxyz ( 1:nx   ,  1:ny   ,  1:nz  )

  fx   (-1:nx+2 , -1:ny+2 , -1:nz+2) = 0.0
  fy   (-1:nx+2 , -1:ny+2 , -1:nz+2) = 0.0
  fz   (-1:nx+2 , -1:ny+2 , -1:nz+2) = 0.0
  fxy  ( 0:nx+1 ,  0:ny+1 ,  0:nz+1) = 0.0
  fxz  ( 0:nx+1 ,  0:ny+1 ,  0:nz+1) = 0.0
  fyz  ( 0:nx+1 ,  0:ny+1 ,  0:nz+1) = 0.0
  fxyz ( 1:nx   ,  1:ny   ,  1:nz  ) = 0.0

  return
end subroutine Interpolate_cubic3DmonoDerv
