!!****if* source/numericalTools/Interpolate/InterpolateMain/Interpolate_cubic3DmonoDerv
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

  integer :: i, j, k

  real    :: Bxy, Bxz, Byz, Bxyz
  real    :: d21, d31, d41, d52, d53, d62, d64, d73, d74, d85, d86, d87
  real    :: dm, dp
  real    :: f000, fp00, f0p0, f00p, fpp0, fp0p, f0pp, fppp
  real    :: fx000, fx0p0, fx00p, fx0m0, fx00m, &
             fy000, fyp00, fy00p, fym00, fy00m, &
             fz000, fzp00, fz0p0, fzm00, fz0m0
  real    :: fxyEstimate, fxzEstimate, fyzEstimate, fxyzEstimate
  real    :: fxy000, fxy00p, fxy00m, fxz000, fxz0p0, fxz0m0, fyz000, fyzp00, fyzm00

  real, parameter :: zero   = 0.0
  real, parameter :: sixth  = 1.0 / 6.0
  real, parameter :: fourth = 0.25
  real, parameter :: half   = 0.5
!
!
!     ...Calculate initial 1st order x,y and z derivatives for monotonic border.
!
!        For each cube on the grid, the following 1st order derivatives are
!        adjusted along the indicated borders:
!
!
!
!                            ----1-----
!                          /|         /|        1 --> fx
!                         2 |        2 |
!                        /  3       /  3        2 --> fy
!                        ----1-----    |
!                       |   |      |   |        3 --> fz
!                       |    ----1-|---
!                       3  /       3  /
!                       | 2        | 2
!                       |/         |/
!                        ----1-----
!
!
!
  do k = -1, nz+2
     do j = -1, ny+2
        do i = -1, nx+2

           f000 = f (i,j,k)

           dm = f000 - f (i-1,j,k)
           dp = f (i+1,j,k) - f000

           if (dm < zero .and. dp < zero) then
               fx (i,j,k) = max (dm,dp)
           else if (dm > zero .and. dp > zero) then
               fx (i,j,k) = min (dm,dp)
           else
               fx (i,j,k) = zero
           end if

           dm = f000 - f (i,j-1,k)
           dp = f (i,j+1,k) - f000

           if (dm < zero .and. dp < zero) then
               fy (i,j,k) = max (dm,dp)
           else if (dm > zero .and. dp > zero) then
               fy (i,j,k) = min (dm,dp)
           else
               fy (i,j,k) = zero
           end if

           dm = f000 - f (i,j,k-1)
           dp = f (i,j,k+1) - f000

           if (dm < zero .and. dp < zero) then
               fz (i,j,k) = max (dm,dp)
           else if (dm > zero .and. dp > zero) then
               fz (i,j,k) = min (dm,dp)
           else
               fz (i,j,k) = zero
           end if

        end do
     end do
  end do
!
!
!     ...Loop over all cubes and set appropriate 1st order derivatives to zero,
!        if cube is warped in x,y and/or z direction.
!
!
  do k = -1, nz+1
     do j = -1, ny+1
        do i = -1, nx+1

           f000 = f (i  ,j  ,k  )
           fp00 = f (i+1,j  ,k  )
           f0p0 = f (i  ,j+1,k  )
           fpp0 = f (i+1,j+1,k  )
           f00p = f (i  ,j  ,k+1)
           fp0p = f (i+1,j  ,k+1)
           f0pp = f (i  ,j+1,k+1)
           fppp = f (i+1,j+1,k+1)

           d21 = fp00 - f000
           d53 = fpp0 - f0p0
           d64 = fp0p - f00p
           d87 = fppp - f0pp

           if (d21 * d53 < zero .or. &
               d21 * d64 < zero .or. &
               d21 * d87 < zero .or. &
               d53 * d64 < zero .or. &
               d53 * d87 < zero .or. &
               d64 * d87 < zero      ) then

               fx (i  ,j  ,k  ) = zero
               fx (i+1,j  ,k  ) = zero
               fx (i  ,j+1,k  ) = zero
               fx (i+1,j+1,k  ) = zero
               fx (i  ,j  ,k+1) = zero
               fx (i+1,j  ,k+1) = zero
               fx (i  ,j+1,k+1) = zero
               fx (i+1,j+1,k+1) = zero
           end if

           d31 = f0p0 - f000
           d52 = fpp0 - fp00
           d74 = f0pp - f00p
           d86 = fppp - fp0p

           if (d31 * d52 < zero .or. &
               d31 * d74 < zero .or. &
               d31 * d86 < zero .or. &
               d52 * d74 < zero .or. &
               d52 * d86 < zero .or. &
               d74 * d86 < zero      ) then

               fy (i  ,j  ,k  ) = zero
               fy (i+1,j  ,k  ) = zero
               fy (i  ,j+1,k  ) = zero
               fy (i+1,j+1,k  ) = zero
               fy (i  ,j  ,k+1) = zero
               fy (i+1,j  ,k+1) = zero
               fy (i  ,j+1,k+1) = zero
               fy (i+1,j+1,k+1) = zero
           end if

           d41 = f00p - f000
           d62 = fp0p - fp00
           d73 = f0pp - f0p0
           d85 = fppp - fpp0

           if (d41 * d62 < zero .or. &
               d41 * d73 < zero .or. &
               d41 * d85 < zero .or. &
               d62 * d73 < zero .or. &
               d62 * d85 < zero .or. &
               d73 * d85 < zero      ) then

               fz (i  ,j  ,k  ) = zero
               fz (i+1,j  ,k  ) = zero
               fz (i  ,j+1,k  ) = zero
               fz (i+1,j+1,k  ) = zero
               fz (i  ,j  ,k+1) = zero
               fz (i+1,j  ,k+1) = zero
               fz (i  ,j+1,k+1) = zero
               fz (i+1,j+1,k+1) = zero
           end if

        end do
     end do
  end do
!
!
!     ...Apply upsweep in x,y and z directions on 1st order derivatives.
!        From grid point 1 (the loop i,j,k basis) do:
!
!
!                            ----------
!                          /|         /|        at 2 --> damp fy and fz
!                         / |        / |
!                        /  |       /  |        at 3 --> damp fx and fz
!                       4----------    |
!                       |   |      |   |        at 4 --> damp fx and fy
!                       |   3------|---
!                       |  /       |  /
!                       | /        | /
!                       |/         |/
!                       1----------2
!
!        Note, that except for the boundary points, each fx,fy,fz inside the
!        grid is addressed (damped) twice (from loop points *) :
!
!                      fx      *----------fy      *----------fz
!                      /|                 |                 /
!                     / |                 |                /
!                    /  |                 |               /
!                   *   |                 |              *
!                       |                 |
!                       *                 *
!
!
  do k = -1, nz+1
     do j = -1, ny+1
        do i = -1, nx+1

           d21 = half * abs (f (i+1,j  ,k  ) - f (i,j,k))
           d31 = half * abs (f (i  ,j+1,k  ) - f (i,j,k))
           d41 = half * abs (f (i  ,j  ,k+1) - f (i,j,k))

           fx000 = abs (fx (i  ,j  ,k  ))
           fx0p0 = abs (fx (i  ,j+1,k  ))
           fx00p = abs (fx (i  ,j  ,k+1))
           fy000 = abs (fy (i  ,j  ,k  ))
           fyp00 = abs (fy (i+1,j  ,k  ))
           fy00p = abs (fy (i  ,j  ,k+1))
           fz000 = abs (fz (i  ,j  ,k  ))
           fzp00 = abs (fz (i+1,j  ,k  ))
           fz0p0 = abs (fz (i  ,j+1,k  ))

           fx0p0 = min (fx0p0, fx000 + d31)
           fx00p = min (fx00p, fx000 + d41)
           fyp00 = min (fyp00, fy000 + d21)
           fy00p = min (fy00p, fy000 + d41)
           fzp00 = min (fzp00, fz000 + d21)
           fz0p0 = min (fz0p0, fz000 + d31)

           fx (i  ,j+1,k  ) = sign (fx0p0, fx (i  ,j+1,k  ))
           fx (i  ,j  ,k+1) = sign (fx00p, fx (i  ,j  ,k+1))
           fy (i+1,j  ,k  ) = sign (fyp00, fy (i+1,j  ,k  ))
           fy (i  ,j  ,k+1) = sign (fy00p, fy (i  ,j  ,k+1))
           fz (i+1,j  ,k  ) = sign (fzp00, fz (i+1,j  ,k  ))
           fz (i  ,j+1,k  ) = sign (fz0p0, fz (i  ,j+1,k  ))

        end do
     end do
  end do
!
!
!     ...Apply downsweep in x,y and z directions on 1st order derivatives.
!        From grid point 8 (the loop i,j,k basis) do:
!
!
!                           7----------8
!                          /|         /|        at 5 --> damp fx and fy
!                         / |        / |
!                        /  |       /  |        at 6 --> damp fx and fz
!                        ----------6   |
!                       |   |      |   |        at 7 --> damp fy and fz
!                       |    ------|---5
!                       |  /       |  /
!                       | /        | /
!                       |/         |/
!                        ----------
!
!        Again, except for the boundary points, each fx,fy,fz inside the
!        grid is addressed (damped) twice (from loop points *) :
!
!                       *          *
!                       |          |
!                       |   *      |                    *
!                       |  /       |                   /
!                       | /        |                  /
!                       |/         |                 /
!                       fx        fy----------*    fz----------*
!
!
  do k = nz+2, 0, -1
     do j = ny+2, 0, -1
        do i = nx+2, 0, -1

           d87 = half * abs (f (i,j,k) - f (i-1,j  ,k  ))
           d86 = half * abs (f (i,j,k) - f (i  ,j-1,k  ))
           d85 = half * abs (f (i,j,k) - f (i  ,j  ,k-1))

           fx000 = abs (fx (i  ,j  ,k  ))
           fx0m0 = abs (fx (i  ,j-1,k  ))
           fx00m = abs (fx (i  ,j  ,k-1))
           fy000 = abs (fy (i  ,j  ,k  ))
           fym00 = abs (fy (i-1,j  ,k  ))
           fy00m = abs (fy (i  ,j  ,k-1))
           fz000 = abs (fz (i  ,j  ,k  ))
           fzm00 = abs (fz (i-1,j  ,k  ))
           fz0m0 = abs (fz (i  ,j-1,k  ))

           fx0m0 = min (fx0m0, fx000 + d86)
           fx00m = min (fx00m, fx000 + d85)
           fym00 = min (fym00, fy000 + d87)
           fy00m = min (fy00m, fy000 + d85)
           fzm00 = min (fzm00, fz000 + d87)
           fz0m0 = min (fz0m0, fz000 + d86)

           fx (i  ,j-1,k  ) = sign (fx0m0, fx (i  ,j-1,k  ))
           fx (i  ,j  ,k-1) = sign (fx00m, fx (i  ,j  ,k-1))
           fy (i-1,j  ,k  ) = sign (fym00, fy (i-1,j  ,k  ))
           fy (i  ,j  ,k-1) = sign (fy00m, fy (i  ,j  ,k-1))
           fz (i-1,j  ,k  ) = sign (fzm00, fz (i-1,j  ,k  ))
           fz (i  ,j-1,k  ) = sign (fz0m0, fz (i  ,j-1,k  ))

        end do
     end do
  end do
!
!
!     ...Calculate the 2nd mixed order derivatives at each grid point.
!
!
  do k = 0,nz+1
     do j = 0,ny+1
        do i = 0,nx+1

           fxyEstimate = fourth * (   fx (i  ,j+1,k  ) &
                                    - fx (i  ,j-1,k  ) &
                                    + fy (i+1,j  ,k  ) &
                                    - fy (i-1,j  ,k  ) )

           fxzEstimate = fourth * (   fx (i  ,j  ,k+1) &
                                    - fx (i  ,j  ,k-1) &
                                    + fz (i+1,j  ,k  ) &
                                    - fz (i-1,j  ,k  ) )

           fyzEstimate = fourth * (   fy (i  ,j  ,k+1) &
                                    - fy (i  ,j  ,k-1) &
                                    + fz (i  ,j+1,k  ) &
                                    - fz (i  ,j-1,k  ) )

           Bxy = min (abs (fx (i,j,k)) , abs (fy (i,j,k)))  ! upper limit for fxy
           Bxz = min (abs (fx (i,j,k)) , abs (fz (i,j,k)))  ! upper limit for fxz
           Byz = min (abs (fy (i,j,k)) , abs (fz (i,j,k)))  ! upper limit for fyz

           fxy (i,j,k) = min (Bxy , max (-Bxy, fxyEstimate))
           fxz (i,j,k) = min (Bxz , max (-Bxz, fxzEstimate))
           fyz (i,j,k) = min (Byz , max (-Byz, fyzEstimate))

        end do
     end do
  end do
!
!
!     ...Apply sign-correction algorithm to the 2nd mixed order derivatives at each
!        relevant grid direction:
!
!                    fxy -> sign-correction along z-direction
!                    fxz -> sign-correction along y-direction
!                    fyz -> sign-correction along x-direction
!
!
  do k = 0,nz
     do j = 0,ny
        do i = 0,nx

           fxy000 = fxy (i,j,k  )
           fxy00p = fxy (i,j,k+1)

           if (fxy000 * fxy00p < zero) then
               if (abs (fxy000) < abs (fxy00p)) then
                   fxy (i,j,k)   = zero
               else if (abs (fxy000) > abs (fxy00p)) then
                   fxy (i,j,k+1) = zero
               else
                   fxy (i,j,k)   = zero
                   fxy (i,j,k+1) = zero
               end if
           end if

           fxz000 = fxz (i,j  ,k)
           fxz0p0 = fxz (i,j+1,k)

           if (fxz000 * fxz0p0 < zero) then
               if (abs (fxz000) < abs (fxz0p0)) then
                   fxz (i,j,k)   = zero
               else if (abs (fxz000) > abs (fxz0p0)) then
                   fxz (i,j+1,k) = zero
               else
                   fxz (i,j,k)   = zero
                   fxz (i,j+1,k) = zero
               end if
           end if

           fyz000 = fyz (i  ,j,k)
           fyzp00 = fyz (i+1,j,k)

           if (fyz000 * fyzp00 < zero) then
               if (abs (fyz000) < abs (fyzp00)) then
                   fyz (i,j,k)   = zero
               else if (abs (fyz000) > abs (fyzp00)) then
                   fyz (i+1,j,k) = zero
               else
                   fyz (i,j,k)   = zero
                   fyz (i+1,j,k) = zero
               end if
           end if


        end do
     end do
  end do
!
!
!     ...Apply upsweep in x,y and z directions on 2nd order derivatives.
!        From grid point 1 (the loop i,j,k basis) do:
!
!
!                            ----------
!                          /|         /|        at 2 --> damp fyz
!                         / |        / |
!                        /  |       /  |        at 3 --> damp fxz
!                       4----------    |
!                       |   |      |   |        at 4 --> damp fxy
!                       |   3------|---
!                       |  /       |  /
!                       | /        | /
!                       |/         |/
!                       1----------2
!
!
  do k = 0, nz
     do j = 0, ny
        do i = 0, nx

           d21 = half * abs (f (i+1,j  ,k  ) - f (i,j,k))
           d31 = half * abs (f (i  ,j+1,k  ) - f (i,j,k))
           d41 = half * abs (f (i  ,j  ,k+1) - f (i,j,k))

           fxy000 = abs (fxy (i  ,j  ,k  ))
           fxy00p = abs (fxy (i  ,j  ,k+1))
           fxz000 = abs (fxz (i  ,j  ,k  ))
           fxz0p0 = abs (fxz (i  ,j+1,k  ))
           fyz000 = abs (fyz (i  ,j  ,k  ))
           fyzp00 = abs (fyz (i+1,j  ,k  ))

           fxy00p = min (fxy00p, fxy000 + d41)
           fxz0p0 = min (fxz0p0, fxz000 + d31)
           fyzp00 = min (fyzp00, fyz000 + d21)

           fxy (i  ,j  ,k+1) = sign (fxy00p, fxy (i  ,j  ,k+1))
           fxz (i  ,j+1,k  ) = sign (fxz0p0, fxz (i  ,j+1,k  ))
           fyz (i+1,j  ,k  ) = sign (fyzp00, fyz (i+1,j  ,k  ))

        end do
     end do
  end do
!
!
!     ...Apply downsweep in x,y and z directions on 2nd order derivatives.
!        From grid point 8 (the loop i,j,k basis) do:
!
!
!                           7----------8
!                          /|         /|        at 5 --> damp fxy
!                         / |        / |
!                        /  |       /  |        at 6 --> damp fxz
!                        ----------6   |
!                       |   |      |   |        at 7 --> damp fyz
!                       |    ------|---5
!                       |  /       |  /
!                       | /        | /
!                       |/         |/
!                        ----------
!
!
  do k = nz+1, 1, -1
     do j = ny+1, 1, -1
        do i = nx+1, 1, -1

           d87 = half * abs (f (i,j,k) - f (i-1,j  ,k  ))
           d86 = half * abs (f (i,j,k) - f (i  ,j-1,k  ))
           d85 = half * abs (f (i,j,k) - f (i  ,j  ,k-1))

           fxy000 = abs (fxy (i  ,j  ,k  ))
           fxy00m = abs (fxy (i  ,j  ,k-1))
           fxz000 = abs (fxz (i  ,j  ,k  ))
           fxz0m0 = abs (fxz (i  ,j-1,k  ))
           fyz000 = abs (fyz (i  ,j  ,k  ))
           fyzm00 = abs (fyz (i-1,j  ,k  ))

           fxy00m = min (fxy00m, fxy000 + d85)
           fxz0m0 = min (fxz0m0, fxz000 + d86)
           fyzm00 = min (fyzm00, fyz000 + d87)

           fxy (i  ,j  ,k-1) = sign (fxy00m, fxy (i  ,j  ,k-1))
           fxz (i  ,j-1,k  ) = sign (fxz0m0, fxz (i  ,j-1,k  ))
           fyz (i-1,j  ,k  ) = sign (fyzm00, fyz (i-1,j  ,k  ))

        end do
     end do
  end do
!
!
!     ...Calculate the 3rd mixed order derivatives at each grid point.
!
!
  do k = 1,nz
     do j = 1,ny
        do i = 1,nx

           fxyzEstimate = sixth * (   fxy (i  ,j  ,k+1) &
                                    - fxy (i  ,j  ,k-1) &
                                    + fxz (i  ,j+1,k  ) &
                                    - fxz (i  ,j-1,k  ) &
                                    + fyz (i+1,j  ,k  ) &
                                    - fyz (i-1,j  ,k  ) )

           Bxyz = min (abs (fx (i,j,k)) , abs (fy (i,j,k)) , abs (fz (i,j,k)))  ! upper limit for fxyz

           fxyz (i,j,k) = min (Bxyz , max (-Bxyz, fxyzEstimate))

        end do
     end do
  end do
!
!
!     ...Ready!
!
!
  return
end subroutine Interpolate_cubic3DmonoDerv
