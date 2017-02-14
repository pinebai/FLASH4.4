!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/Gravity_bhEwaldPotV42
!!
!! NAME
!!
!!  Gravity_bhEwald
!!
!!
!! SYNOPSIS
!!
!!   real field = Gravity_bhEwaldPotV42(
!!                           real(in)    :: x,
!!                           real(in)    :: y,
!!                           real(in)    :: z
!!        )
!!
!! DESCRIPTION
!!
!!   Interpolates in the Ewald field and returns its value for the point x,y,z.
!!
!! ARGUMENTS
!!
!!  x   : x-coordinate of the point where the Ewald field is determined
!!  y   : y-coordinate of the point where the Ewald field is determined
!!  z   : z-coordinate of the point where the Ewald field is determined
!!
!! RESULT
!!
!!  Value of the Ewald field in a given point.
!!
!! NOTES
!!
!!***



function Gravity_bhEwaldPotV42(x, y, z)

  use Gravity_data, ONLY: grv_bhEwaldFieldNxV42, grv_bhEwaldFieldNyV42, grv_bhEwaldFieldNzV42, &
      & grv_bhTreeEwaldPotV42, grv_bhEwaldLMaxV42, grv_bhDxI, grv_bhDyIV42, grv_bhDzIV42, &
      & grv_bhMinEFSizeV42, grv_bhEwaldNRefV42, grv_bhDirectionQuadV42
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  real,intent(in) :: x, y, z
  real :: Gravity_bhEwaldPotV42
  real :: y0, y1, y2
  real :: dx, dy, dz, p, q, r
  real :: dist_max
  real :: wf(1:8)
  integer :: i, j, k, m
  integer :: ip0, ip1, ip2, jp0, jp1, jp2, kp0, kp1, kp2
  integer :: ewald_reflevel, c, ncount

  ! choose which level in Ewald field refinement to interpolate in
  dist_max = max(x,y,z)
  ! interpolation deeper than 0.5 in the hightest level means that 
  ! higher value of grv_bhEwaldNRefV42 must be choosen in flash.par
  if (dist_max <= grv_bhMinEFSizeV42) then
    ewald_reflevel = grv_bhEwaldNRefV42-1
    if (dist_max < grv_bhMinEFSizeV42) then
      write(*,*) 'Interpolation is too deep in the highest level ' &
      &         ,'ratio is ',dist_max/grv_bhMinEFSizeV42
      call Driver_abortFlash("FATAL: Interpolation is unreliable")
    end if
  else
    ewald_reflevel = -1
    do while (dist_max < grv_bhEwaldLMaxV42)
      dist_max = dist_max*2
      ewald_reflevel = ewald_reflevel + 1
    end do
  end if

  ! find the nearest smaller indeces in a given field
  c = 2**ewald_reflevel
  i = floor(x*grv_bhDxI*c)
  j = floor(y*grv_bhDyIV42*c)
  k = floor(z*grv_bhDzIV42*c)

  if ((i < 0) .or. (i >= grv_bhEwaldFieldNxV42)) then
    write(*,*) 'Ewald: ',i,j,k,ewald_reflevel, 1,c,grv_bhDxI,grv_bhDyIV42 &
    & ,grv_bhDzIV42,dist_max,grv_bhMinEFSizeV42, x, y, z, grv_bhEwaldLMaxV42
    call Driver_abortFlash("FATAL: Ewald field i out of limits")
  endif
  if ((j < 0) .or. (j >= grv_bhEwaldFieldNyV42)) then
    write(*,*) 'Ewald: ',i,j,k,ewald_reflevel, 1,c,grv_bhDxI,grv_bhDyIV42 &
    & ,grv_bhDzIV42,dist_max,grv_bhMinEFSizeV42
    call Driver_abortFlash("FATAL: Ewald field j out of limits")
  endif
  if ((k < 0) .or. (k >= grv_bhEwaldFieldNzV42)) then
    write(*,*) 'Ewald: ',i,j,k,ewald_reflevel, 1,c,grv_bhDxI,grv_bhDyIV42 &
    & ,grv_bhDzIV42,dist_max,grv_bhMinEFSizeV42
    call Driver_abortFlash("FATAL: Ewald field k out of limits")
  endif

  p = x*grv_bhDxI*c - i
  q = y*grv_bhDyIV42*c - j
  r = z*grv_bhDzIV42*c - k

  ip1 = i+1
  jp1 = j+1
  kp1 = k+1

  Select case (grv_bhDirectionQuadV42)
! linear interpolation
  Case (0)

    wf(1)=(1.0 - p)*(1.0 - q)*(1.0 - r)
    wf(2)=(1.0 - p)*(1.0 - q)*(      r)
    wf(3)=(1.0 - p)*(      q)*(1.0 - r)
    wf(4)=(1.0 - p)*(      q)*(      r)
    wf(5)=(      p)*(1.0 - q)*(1.0 - r)
    wf(6)=(      p)*(1.0 - q)*(      r)
    wf(7)=(      p)*(      q)*(1.0 - r)
    wf(8)=(      p)*(      q)*(      r)

    Gravity_bhEwaldPotV42 = ( &
      & wf(1) * grv_bhTreeEwaldPotV42(ewald_reflevel, i  , j  , k  ) + &
      & wf(2) * grv_bhTreeEwaldPotV42(ewald_reflevel, i  , j  , kp1) + &
      & wf(3) * grv_bhTreeEwaldPotV42(ewald_reflevel, i  , jp1, k  ) + &
      & wf(4) * grv_bhTreeEwaldPotV42(ewald_reflevel, i  , jp1, kp1) + &
      & wf(5) * grv_bhTreeEwaldPotV42(ewald_reflevel, ip1, j  , k  ) + &
      & wf(6) * grv_bhTreeEwaldPotV42(ewald_reflevel, ip1, j  , kp1) + &
      & wf(7) * grv_bhTreeEwaldPotV42(ewald_reflevel, ip1, jp1, k  ) + &
      & wf(8) * grv_bhTreeEwaldPotV42(ewald_reflevel, ip1, jp1, kp1)   &
      & )

! quadratic interpolation in direction x
    Case (1)
  ! note: Ewald_field stretches from -1 to grv_bhEwaldFieldNxV42
  if (p.ge.0.5) then
    ip0 = i
    ip1 = i+1
    ip2 = i+2
    p = p-1.0
  else
    ip0 = i-1
    ip1 = i
    ip2 = i+1
  endif

  ! evaluate the Ewald field in three subsequent slices parallel to plane x=0  
  y0 = ( &
    & (1.0 - q)*(1.0 - r)*grv_bhTreeEwaldPotV42(ewald_reflevel, ip0, j  , k    ) + &
    & (1.0 - q)*(      r)*grv_bhTreeEwaldPotV42(ewald_reflevel, ip0, j  , kp1  ) + &
    & (      q)*(1.0 - r)*grv_bhTreeEwaldPotV42(ewald_reflevel, ip0, jp1, k    ) + &
    & (      q)*(      r)*grv_bhTreeEwaldPotV42(ewald_reflevel, ip0, jp1, kp1  )   &
    & )

  y1 = ( &
    & (1.0 - q)*(1.0 - r)*grv_bhTreeEwaldPotV42(ewald_reflevel, ip1, j  , k    ) + &
    & (1.0 - q)*(      r)*grv_bhTreeEwaldPotV42(ewald_reflevel, ip1, j  , kp1  ) + &
    & (      q)*(1.0 - r)*grv_bhTreeEwaldPotV42(ewald_reflevel, ip1, jp1, k    ) + &
    & (      q)*(      r)*grv_bhTreeEwaldPotV42(ewald_reflevel, ip1, jp1, kp1  )   &
    & )

  y2 = ( &
    & (1.0 - q)*(1.0 - r)*grv_bhTreeEwaldPotV42(ewald_reflevel, ip2, j  , k    ) + &
    & (1.0 - q)*(      r)*grv_bhTreeEwaldPotV42(ewald_reflevel, ip2, j  , kp1  ) + &
    & (      q)*(1.0 - r)*grv_bhTreeEwaldPotV42(ewald_reflevel, ip2, jp1, k    ) + &
    & (      q)*(      r)*grv_bhTreeEwaldPotV42(ewald_reflevel, ip2, jp1, kp1  )   &
    & )

  Gravity_bhEwaldPotV42 = y1+0.5*(p*p*(y0+y2-2.0*y1)+p*(y2-y0))

! quadratic interpolation in direction y
    Case (2)
  ! note: Ewald_field stretches from -1 to grv_bhEwaldFieldNxV42
  if (q.ge.0.5) then
    jp0 = j
    jp1 = j+1
    jp2 = j+2
    q = q-1.0
  else
    jp0 = j-1
    jp1 = j
    jp2 = j+1
  endif

  ! evaluate the Ewald field in three subsequent slices parellel to plane z=0  
  y0 = ( &
    & (1.0 - p)*(1.0 - r)*grv_bhTreeEwaldPotV42(ewald_reflevel, i  , jp0, k    ) + &
    & (1.0 - p)*(      r)*grv_bhTreeEwaldPotV42(ewald_reflevel, i  , jp0, kp1  ) + &
    & (      p)*(1.0 - r)*grv_bhTreeEwaldPotV42(ewald_reflevel, ip1, jp0, k    ) + &
    & (      p)*(      r)*grv_bhTreeEwaldPotV42(ewald_reflevel, ip1, jp0, kp1  )   &
    & )

  y1 = ( &
    & (1.0 - p)*(1.0 - r)*grv_bhTreeEwaldPotV42(ewald_reflevel, i  , jp1, k    ) + &
    & (1.0 - p)*(      r)*grv_bhTreeEwaldPotV42(ewald_reflevel, i  , jp1, kp1  ) + &
    & (      p)*(1.0 - r)*grv_bhTreeEwaldPotV42(ewald_reflevel, ip1, jp1, k    ) + &
    & (      p)*(      r)*grv_bhTreeEwaldPotV42(ewald_reflevel, ip1, jp1, kp1  )   &
    & )

  y2 = ( &
    & (1.0 - p)*(1.0 - r)*grv_bhTreeEwaldPotV42(ewald_reflevel, i  , jp2, k    ) + &
    & (1.0 - p)*(      r)*grv_bhTreeEwaldPotV42(ewald_reflevel, i  , jp2, kp1  ) + &
    & (      p)*(1.0 - r)*grv_bhTreeEwaldPotV42(ewald_reflevel, ip1, jp2, k    ) + &
    & (      p)*(      r)*grv_bhTreeEwaldPotV42(ewald_reflevel, ip1, jp2, kp1  )   &
    & )

  Gravity_bhEwaldPotV42 = y1+0.5*(q*q*(y0+y2-2.0*y1)+q*(y2-y0))

! quadratic interpolation in direction z
    Case (3)
  ! note: Ewald_field stretches from -1 to grv_bhEwaldFieldNxV42
  if (r.ge.0.5) then
    kp0 = k
    kp1 = k+1
    kp2 = k+2
    r = r-1.0
  else
    kp0 = k-1
    kp1 = k
    kp2 = k+1
  endif

  ! evaluate the Ewald field in three subsequent slices parellel to plane z=0  
  y0 = ( &
    & (1.0 - p)*(1.0 - q)*grv_bhTreeEwaldPotV42(ewald_reflevel, i  , j  , kp0  ) + &
    & (1.0 - p)*(      q)*grv_bhTreeEwaldPotV42(ewald_reflevel, i  , jp1, kp0  ) + &
    & (      p)*(1.0 - q)*grv_bhTreeEwaldPotV42(ewald_reflevel, ip1, j  , kp0  ) + &
    & (      p)*(      q)*grv_bhTreeEwaldPotV42(ewald_reflevel, ip1, jp1, kp0  )   &
    & )

  y1 = ( &
    & (1.0 - p)*(1.0 - q)*grv_bhTreeEwaldPotV42(ewald_reflevel, i  , j  , kp1  ) + &
    & (1.0 - p)*(      q)*grv_bhTreeEwaldPotV42(ewald_reflevel, i  , jp1, kp1  ) + &
    & (      p)*(1.0 - q)*grv_bhTreeEwaldPotV42(ewald_reflevel, ip1, j  , kp1  ) + &
    & (      p)*(      q)*grv_bhTreeEwaldPotV42(ewald_reflevel, ip1, jp1, kp1  )   &
    & )

  y2 = ( &
    & (1.0 - p)*(1.0 - q)*grv_bhTreeEwaldPotV42(ewald_reflevel, i  , j  , kp2  ) + &
    & (1.0 - p)*(      q)*grv_bhTreeEwaldPotV42(ewald_reflevel, i  , jp1, kp2  ) + &
    & (      p)*(1.0 - q)*grv_bhTreeEwaldPotV42(ewald_reflevel, ip1, j  , kp2  ) + &
    & (      p)*(      q)*grv_bhTreeEwaldPotV42(ewald_reflevel, ip1, jp1, kp2  )   &
    & )

  Gravity_bhEwaldPotV42 = y1+0.5*(r*r*(y0+y2-2.0*y1)+r*(y2-y0))


  End select

  return
end

