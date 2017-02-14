!!****if* source/physics/Hydro/HydroMain/split/RHD/hy_rhd_hlle
!!
!! NAME
!!
!!   hy_rhd_hlle
!!
!! SYNOPSIS
!!   hy_rhd_hlle( real   (IN) :: Vc(:),
!!                real   (IN) :: Vm(:),
!!                real   (IN) :: Vp(:),
!!                real  (OUT) :: Flux(:,:),
!!                real   (IN) :: x(:),
!!                real   (IN) :: dx(:),
!!                real   (IN) :: dt,
!!                real  (OUT) :: speed(:),
!!                real  (OUT) :: vint(:),
!!                integer(IN) :: ibeg,
!!                integer(IN) :: iend,
!!                integer(IN) :: n,
!!                integer(IN) :: dir)
!!
!!
!! DESCRIPTION
!!
!!   Solve Riemann problem using a two-wave HLLE solver. 
!!   Return numerical fluxes.
!!
!!
!! ARGUMENTS
!! 
!!   Vc     -      Array of centred values
!!   Vm,Vp  -      Arrays of left and right interpolated variables
!!   Flux   -      Array of numerical fluxes
!!   speed  -      On ouput, it contains the fastest characteristic 
!!                 speed computed at the solution of the Riemann problem
!!   x      -      locations of cell centroids.    
!!   dx     -      mesh spacing 
!!   dt     -      time step
!!   vint   -      Array of interface velocities passed to the calling
!!                 routine that advances species in Lagrangian fashion
!!   ibeg   -      initial point
!!   iend   -      final point
!!   n      -      Size of arrays in the sweep direction
!!   dir    -      Sweep direction
!!
!!***
subroutine hy_rhd_hlle(Vc, Vm, Vp, Flux, x, dx, dt, speed, vint, ibeg, iend, n, dir)

  use hy_rhd_interface, ONLY: hy_rhd_flux,      &
                              hy_rhd_enthalpy,  &
                              hy_rhd_maxChSpeed,&
                              hy_rhd_soundSpeed2

  implicit none

#include "Flash.h"
#include "constants.h"
#include "RHD.h"

  integer, INTENT(in) :: n, dir, ibeg, iend
  real, DIMENSION(NUNK_VARS,n), INTENT(in)  :: Vm, Vp, Vc
  real, DIMENSION(NFLUXES,n), INTENT(out) :: Flux
  real, DIMENSION(n), INTENT(in)  :: x, dx
  real, DIMENSION(n), INTENT(out) :: speed, vint
  real, INTENT(in)                :: dt

  integer                  ::  i, nv
  real, DIMENSION(NUNK_VARS,n)  ::  Up, Um, Vs
  real, DIMENSION(NFLUXES,n)  ::  Fm, Fp
  real, DIMENSION(n)       ::  hs 
  real, DIMENSION(NUNK_VARS)    ::  Vstar, qR, qL

  real, DIMENSION(n)       :: hp, hm, am2, ap2  
  real :: cs2r, cs2l, vxr, vxl, d2r, d2l
  real :: vt2r, vt2l, vel2r, vel2l, scrhr, scrhl, bmax, bmin
 
  integer  :: VELN_VAR, VELT1_VAR, VELT2_VAR,&
              VELN_FLUX,VELT1_FLUX,VELT2_FLUX


  select case(dir)
  case (SWEEP_X)
     VELN_VAR = VELX_VAR
     VELT1_VAR= VELY_VAR
     VELT2_VAR= VELZ_VAR

     VELN_FLUX = XMOM_FLUX
     VELT1_FLUX= YMOM_FLUX
     VELT2_FLUX= ZMOM_FLUX
  case (SWEEP_Y)
     VELN_VAR = VELY_VAR
     VELT1_VAR= VELZ_VAR
     VELT2_VAR= VELX_VAR

     VELN_FLUX = YMOM_FLUX
     VELT1_FLUX= ZMOM_FLUX
     VELT2_FLUX= XMOM_FLUX
  case (SWEEP_Z)
     VELN_VAR = VELZ_VAR
     VELT1_VAR= VELX_VAR
     VELT2_VAR= VELY_VAR

     VELN_FLUX = ZMOM_FLUX
     VELT1_FLUX= XMOM_FLUX
     VELT2_FLUX= YMOM_FLUX
  end select


  vint(:) = 0.0                 !DEV: Is this correct??? - KW
  flux = 0.e0

  call hy_rhd_flux(Vm, Fm, Um, ibeg, iend, n, dir)
  call hy_rhd_flux(Vp, Fp, Up, ibeg - 1, iend - 1, n, dir)

  call hy_rhd_enthalpy (Vm, hm, ibeg, iend, n)
  call hy_rhd_enthalpy (Vp, hp, ibeg - 1, iend - 1, n)

  call hy_rhd_soundSpeed2(Vm, hm, am2, ibeg, iend, n)
  call hy_rhd_soundSpeed2(Vp, hp, ap2, ibeg - 1, iend - 1, n)

  do i = ibeg, iend

    qR   = Vm(:,i)
    qL   = Vp(:,i - 1)
    cs2r = am2(i)
    cs2l = ap2(i - 1)     

    vxr   = qR(VELN_VAR)
    vxl   = qL(VELN_VAR)
    vt2r  = qR(VELT1_VAR)*qR(VELT1_VAR) + qR(VELT2_VAR)*qR(VELT2_VAR)
    vt2l  = qL(VELT1_VAR)*qL(VELT1_VAR) + qL(VELT2_VAR)*qL(VELT2_VAR)
    vel2r = vxr*vxr + vt2r
    vel2l = vxl*vxl + vt2l

    scrhr = sqrt(cs2r*(1.e0 - vxr*vxr - vt2r*cs2r)*(1.e0 - vel2r))
    scrhl = sqrt(cs2l*(1.e0 - vxl*vxl - vt2l*cs2l)*(1.e0 - vel2l))
  
    d2r    = 1.e0 - vel2r*cs2r
    d2l    = 1.e0 - vel2l*cs2l
    
    scrhr = scrhr/d2r
    scrhl = scrhl/d2l

    d2r = (1.e0 - cs2r)/d2r
    d2l = (1.e0 - cs2l)/d2l

! ------------------------------------------------- 
!             Davis Estimate
! ------------------------------------------------- 
   
    bmin = min(vxr*d2r - scrhr, vxl*d2l - scrhl)
    bmin = min(0.e0, bmin)

    bmax = max(vxr*d2r + scrhr, vxl*d2l + scrhl)
    bmax = max(0.e0, bmax)

    speed(i) = max(-bmin,bmax)

    scrhr = 1.0/(bmax - bmin)
    scrhl = bmin*bmax*scrhr
    bmax  = bmax*scrhr
    bmin  = bmin*scrhr

    flux(DENS_FLUX,i) = scrhl*(Um(DENS_VAR,i) - Up(DENS_VAR,i - 1))    &
                     + bmax*Fp(DENS_FLUX,i - 1) - bmin*Fm(DENS_FLUX,i)

    flux(VELN_FLUX,i) = scrhl*(Um(VELN_VAR,i) - Up(VELN_VAR,i - 1))    &
                     + bmax*Fp(VELN_FLUX,i - 1) - bmin*Fm(VELN_FLUX,i)

    flux(VELT1_FLUX,i) = scrhl*(Um(VELT1_VAR,i) - Up(VELT1_VAR,i - 1)) &
                     + bmax*Fp(VELT1_FLUX,i - 1) - bmin*Fm(VELT1_FLUX,i)

    flux(VELT2_FLUX,i) = scrhl*(Um(VELT2_VAR,i) - Up(VELT2_VAR,i - 1)) &
                     + bmax*Fp(VELT2_FLUX,i - 1) - bmin*Fm(VELT2_FLUX,i)

    flux(ENER_FLUX,i) = scrhl*(Um(ENER_VAR,i) - Up(ENER_VAR,i - 1))    &
                     + bmax*Fp(ENER_FLUX,i - 1) - bmin*Fm(ENER_FLUX,i)

  end do
  return

  ! =====================================
  !          LAX FRIEDRICHS
  ! =====================================

  call hy_rhd_flux(Vm, Fm, Um, ibeg, iend, n, dir)
  call hy_rhd_flux(Vp, Fp, Up, ibeg - 1, iend - 1, n, dir)

  do i = ibeg, iend
   Vs(:, i) = 0.5d0*(Vp(:, i - 1) + Vm(:, i))
  end do

  call hy_rhd_enthalpy (Vs, hs, ibeg, iend, n)
  call hy_rhd_maxChSpeed(Vs, hs, speed, ibeg, iend, n, dir)

  do i = ibeg, iend
     flux(DENS_FLUX,i) = 0.5d0*(Fp(DENS_FLUX,i-1) + Fm(DENS_FLUX,i) + &
                      speed(i)*(Up(DENS_VAR, i-1) - Um(DENS_VAR, i)))

     flux(XMOM_FLUX,i) = 0.5d0*(Fp(XMOM_FLUX,i-1) + Fm(XMOM_FLUX,i) + &
                      speed(i)*(Up(VELX_VAR, i-1) - Um(VELX_VAR, i)))

     flux(YMOM_FLUX,i) = 0.5d0*(Fp(YMOM_FLUX,i-1) + Fm(YMOM_FLUX,i) + &
                      speed(i)*(Up(VELY_VAR, i-1) - Um(VELY_VAR, i)))

     flux(ZMOM_FLUX,i) = 0.5d0*(Fp(ZMOM_FLUX,i-1) + Fm(ZMOM_FLUX,i) + &
                      speed(i)*(Up(VELZ_VAR, i-1) - Um(VELZ_VAR, i)))

     flux(ENER_FLUX,i) = 0.5d0*(Fp(ENER_FLUX,i-1) + Fm(ENER_FLUX,i) + &
                      speed(i)*(Up(ENER_VAR, i-1) - Um(ENER_VAR, i)))
  end do

end subroutine hy_rhd_hlle







