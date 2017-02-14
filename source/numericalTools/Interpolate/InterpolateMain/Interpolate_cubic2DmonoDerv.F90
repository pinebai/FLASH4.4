!!****if* source/numericalTools/Interpolate/InterpolateMain/Interpolate_cubic2DmonoDerv
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

  integer :: i, j

  real    :: d1, d2
  real    :: dm0, dp0, d0m, d0p
  real    :: dx0m, dx0p, dym0, dyp0
  real    :: f00, fm0, fp0, f0m, f0p, fpp
  real    :: fx00, fx0m, fx0p, fy00, fym0, fyp0
  real    :: fxyEstimate
  real    :: lowerBound, upperBound
  real    :: s1, s2

  real, parameter :: zero   = 0.0
  real, parameter :: fourth = 0.25
  real, parameter :: one    = 1.0
!
!
!     ...Calculate initial 1st order x and y derivatives for monotonic border.
!
!        For each square on the grid, the following 1st order derivatives are
!        adjusted along the indicated borders:
!
!
!                                      fx
!                                  ----------
!                                 |          |               
!                                 |          |
!                              fy |          | fy
!                                 |          |
!                                 |          |
!                                  ----------
!                                      fx
!
!
  do j = 0, ny+1
     do i = 0, nx+1

        fm0 = f (i-1,j  )
        f00 = f (i  ,j  )
        fp0 = f (i+1,j  )
        f0m = f (i  ,j-1)
        f0p = f (i  ,j+1)

        d1 = f00 - fm0
        d2 = fp0 - f00
        s1 = sign (one,d1)
        s2 = sign (one,d2)

        if (s1 /= s2) then
            fx (i,j) = zero
        else if (s1 < zero) then
            fx (i,j) = max (d1,d2)
        else
            fx (i,j) = min (d1,d2)
        end if

        d1 = f00 - f0m
        d2 = f0p - f00
        s1 = sign (one,d1)
        s2 = sign (one,d2)

        if (s1 /= s2) then
            fy (i,j) = zero
        else if (s1 < zero) then
            fy (i,j) = max (d1,d2)
        else
            fy (i,j) = min (d1,d2)
        end if

     end do
  end do
!
!
!     ...Loop over all squares and set appropriate 1st order derivatives to zero,
!        if square is warped in x and/or y direction.
!
!
  do j = 0, ny
     do i = 0, nx

        f00 = f (i  ,j  )
        fp0 = f (i+1,j  )
        f0p = f (i  ,j+1)
        fpp = f (i+1,j+1)

        d1 = f0p - f00
        d2 = fpp - fp0
        s1 = sign (one,d1)
        s2 = sign (one,d2)

        if (s1 /= s2) then
            fx (i  ,j  ) = zero
            fx (i+1,j  ) = zero
            fx (i  ,j+1) = zero
            fx (i+1,j+1) = zero
        end if

        d1 = fp0 - f00
        d2 = fpp - f0p
        s1 = sign (one,d1)
        s2 = sign (one,d2)

        if (s1 /= s2) then
            fy (i  ,j  ) = zero
            fy (i+1,j  ) = zero
            fy (i  ,j+1) = zero
            fy (i+1,j+1) = zero
        end if

     end do
  end do
!
!
!     ...Apply upsweep in both x and y directions on 1st order derivatives.
!
!
  do j = 0, ny
     do i = 0, nx

        d1 = abs (f (i+1,j  ) - f (i,j))
        d2 = abs (f (i  ,j+1) - f (i,j))

        fy00 = abs (fy (i  ,j  ))
        fyp0 = abs (fy (i+1,j  ))
        fx00 = abs (fx (i  ,j  ))
        fx0p = abs (fx (i  ,j+1))
        fyp0 = min (fyp0, fy00 + d1 + d1)
        fx0p = min (fx0p, fx00 + d2 + d2)

        fy (i+1,j  ) = sign (fyp0, fy (i+1,j  ))
        fx (i  ,j+1) = sign (fx0p, fx (i  ,j+1))

     end do
  end do
!
!
!     ...Apply downsweep in both x and y directions on 1st order derivatives.
!
!
  do j = ny, 0, -1
     do i = nx, 0, -1

        d1 = abs (f (i+1,j  ) - f (i,j))
        d2 = abs (f (i  ,j+1) - f (i,j))

        fy00 = abs (fy (i  ,j  ))
        fyp0 = abs (fy (i+1,j  ))
        fx00 = abs (fx (i  ,j  ))
        fx0p = abs (fx (i  ,j+1))
        fy00 = min (fy00, fyp0 + d1 + d1)
        fx00 = min (fx00, fx0p + d2 + d2)

        fy (i,j) = sign (fy00, fy (i,j))
        fx (i,j) = sign (fx00, fx (i,j))

     end do
  end do
!
!
!     ...Calculate the 2nd mixed order derivatives at each grid point.
!
!
  do j = 1,ny
     do i = 1,nx

        fm0  = f  (i-1,j  )
        f00  = f  (i  ,j  )
        fp0  = f  (i+1,j  )
        f0m  = f  (i  ,j-1)
        f0p  = f  (i  ,j+1)

        fx0m = fx (i  ,j-1)
        fx00 = fx (i  ,j  )
        fx0p = fx (i  ,j+1)
        fym0 = fy (i-1,j  )
        fy00 = fy (i  ,j  )
        fyp0 = fy (i+1,j  )

        dp0  = abs (fp0 - f00)
        dm0  = abs (f00 - fm0)
        d0p  = abs (f0p - f00)
        d0m  = abs (f00 - f0m)

        dx0p = fx0p - fx00
        dx0m = fx00 - fx0m
        dyp0 = fyp0 - fy00
        dym0 = fy00 - fym0

        fx00 = abs (fx00)
        fy00 = abs (fy00)

        lowerBound = max (- fx00,                          &
                          - fy00,                          &
                            dyp0 - dp0 - dp0 - dp0 + fx00, &
                            dym0 - dm0 - dm0 - dm0 + fx00, &
                            dx0p - d0p - d0p - d0p + fy00, &
                            dx0m - d0m - d0m - d0m + fy00  )

        lowerBound = lowerBound + lowerBound + lowerBound

        upperBound = min (  fx00,                          &
                            fy00,                          &
                            dyp0 + dp0 + dp0 + dp0 - fx00, &
                            dym0 + dm0 + dm0 + dm0 - fx00, &
                            dx0p + d0p + d0p + d0p - fy00, &
                            dx0m + d0m + d0m + d0m - fy00  )

        upperBound = upperBound + upperBound + upperBound

        if (lowerBound > upperBound) then
            write (*,'(a)'        ) 'fxy bounds are bad!'
            write (*,'(a,2i2)'    ) 'Occured at i,j position: ',i,j
            write (*,'(a,4ES14.6)') 'Lower / upper bounds are: ', lowerBound, upperBound
            stop
        end if

        fxyEstimate = fourth * (fx0p + fyp0 - fx0m - fym0)

        fxy (i,j) = min (upperBound, max (lowerBound, fxyEstimate))

     end do
  end do
!
!
!     ...Ready!
!
!
  return
end subroutine Interpolate_cubic2DmonoDerv
