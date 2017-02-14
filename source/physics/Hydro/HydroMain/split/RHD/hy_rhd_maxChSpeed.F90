!!****if* source/physics/Hydro/HydroMain/split/RHD/hy_rhd_maxChSpeed
!!
!! NAME
!!
!!   hy_rhd_maxChSpeed
!!
!! SYNOPSIS
!!   hy_rhd_maxChSpeed( real    (IN):: u(:,:),
!!                      real    (IN):: h(:),
!!                      real   (OUT):: speed(:),
!!                      integer (IN):: ibeg,
!!                      integer (IN):: iend,
!!                      integer (IN):: n,
!!                      integer (IN):: dir)
!!
!!
!! DESCRIPTION
!!
!!     Defines the maximum propagation speed for the system of
!!     equations.
!! 
!! ARGUMENTS
!! 
!!  u     -   a 1-D vector of primitive quantities
!!  h     -   a 1-D vector containing the enthalpy
!!  speed -   a 1-D vector containing the maximum characteristic speed
!!  ibeg  -   initial point
!!  iend  -   final point
!!  n     -   number of points
!!  dir   -   flux component; dir = 1,2,3 for the 3 axis
!! 
!!***
subroutine hy_rhd_maxChSpeed (u, h, speed, ibeg, iend, n , dir)  

  use hy_rhd_interface, ONLY: hy_rhd_soundSpeed2

  implicit none

#include "constants.h"
#include "Flash.h"
#include "RHD.h"  

  integer, INTENT(in)                 :: ibeg, iend, n, dir
  real, DIMENSION(NUNK_VARS,n), INTENT(in) :: u
  real, DIMENSION(n), INTENT(in)      :: h
  real, DIMENSION(n), INTENT(out)     :: speed

  integer :: i, ivn
  real, DIMENSION(n)  :: a2
  real                :: sroot, vel2, vx

  ivn = VELX_VAR + dir - IAXIS

  call hy_rhd_soundSpeed2(u, h, a2, ibeg, iend, n)

  do i = ibeg, iend
    vx    = u(ivn,i)
    vel2  = u(VELX_VAR,i)*u(VELX_VAR,i) + &
            u(VELY_VAR,i)*u(VELY_VAR,i) + &
            u(VELZ_VAR,i)*u(VELZ_VAR,i)

    sroot = a2(i)*(1.d0 - vel2)*(1.d0 - vel2*a2(i) - vx*vx*(1.d0 - a2(i)))

    sroot = sqrt(sroot)
    speed(i) = (abs(vx)*(1.d0 - a2(i)) + sroot)/(1.d0 - vel2*a2(i))
  end do
end subroutine hy_rhd_maxChSpeed
