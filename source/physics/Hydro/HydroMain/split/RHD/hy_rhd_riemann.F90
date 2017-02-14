!!****if* source/physics/Hydro/HydroMain/split/RHD/hy_rhd_riemann
!!
!! NAME
!!
!!   hy_rhd_riemann
!!
!! SYNOPSIS
!!  hy_rhd_riemann( real    (IN) :: Vc(:),
!!                  real    (IN) :: Vm(:),
!!                  real    (IN) :: Vp(:),
!!                  real  (OUT) :: Flux(:,:),
!!                  real   (IN) :: x(:),
!!                  real   (IN) :: dx(:),
!!                  real   (IN) :: dt,
!!                  real  (OUT) :: speed(:),
!!                  real  (OUT) :: vint(:),
!!                  integer(IN) :: n,
!!                  integer(IN) :: dir)
!!
!! DESCRIPTION
!!
!!   Given zone edge values Vm and Vp, compute interface fluxes
!!   by solving a Riemann problem.
!!
!!   Reference: Mignone et al, 2005
!!
!!
!! ARGUMENTS
!! 
!!   Vc    -      Array of centred values
!!   Vm,Vp -      Arrays of left and right interpolated variables
!!   Flux  -      Array of RHD fluxes
!!   x     -      locations of cell centroids.    
!!   dx    -      mesh spacing 
!!   dt    -      time step
!!   speed -      On ouput, it contains the fastest characteristic 
!!                 speed computed at the solution of the Riemann problem
!!   vint  -      Array of interface velocities passed to the calling
!!                 routine that advances species in Lagrangian fashion.
!!                 THIS IS CURRENTLY NOT IMPLEMENTED, ALWAYS RETURNS 0.
!!   n     -      Size of arrays in the sweep direction
!!   dir   -      Sweep direction
!!
!!***
subroutine hy_rhd_riemann(Vc, Vm, Vp, Flux, x, dx, dt, speed, vint, n, dir)

  use Hydro_data, ONLY : hy_gamma

  use hy_rhd_interface, ONLY: hy_rhd_flux,       &
                              hy_rhd_enthalpy,   &
                              hy_rhd_maxChSpeed, &
                              hy_rhd_hlle,       &
                              hy_rhd_shock,      &
                              hy_rhd_primitiveToConserve

  implicit none

#include "Flash.h"
#include "constants.h"
#include "RHD.h"

  !! Argument list ----------------------------------------------
  integer, INTENT(in) :: n, dir
  real, DIMENSION(NUNK_VARS,n), INTENT(in)  :: Vm, Vp, Vc
  real, DIMENSION(NFLUXES,n), INTENT(out) :: Flux
  real, DIMENSION(n), INTENT(in)  :: x, dx
  real, DIMENSION(n), INTENT(out) :: speed, vint
  real, INTENT(in)                :: dt
  !! ------------------------------------------------------------

  integer, PARAMETER       ::  MAX_ITER = 20
  integer                  ::  i, iter, nv, nfail
  integer, DIMENSION(n)    ::  failzone
  real, DIMENSION(NUNK_VARS,n)  ::  Vs, Us, Uc
  real, DIMENSION(n)       ::  hp, hm 
  real, DIMENSION(NUNK_VARS)    ::  Vstar, qR, qL
  real  :: duR, gR1, gR, hR, uxR, tauR, uxR1, pR, j2R, zeta_R, vR, vR1
  real  :: duL, gL1, gL, hL, uxL, tauL, uxL1, pL, j2L, zeta_L, vL, vL1
  real  :: p1, dp, a, u1, gamma_r

  integer  :: VELN_VAR,VELT1_VAR,VELT2_VAR
  integer  :: sweepBegin,sweepEnd
 
  !!  call hy_rhd_hlle(Vc, Vm, Vp, Flux, x, dx, dt, speed, vint, 5, n-3, n, dir)
  !!  return

  vint(:) = 0.0

  gamma_r = hy_gamma/(hy_gamma - 1.d0)

  !! define normal and tangent indexes 
  select case(dir)
  case (SWEEP_X)
     VELN_VAR = VELX_VAR
     VELT1_VAR= VELY_VAR
     VELT2_VAR= VELZ_VAR
  case (SWEEP_Y)
     VELN_VAR = VELY_VAR
     VELT1_VAR= VELZ_VAR
     VELT2_VAR= VELX_VAR
  case (SWEEP_Z)
     VELN_VAR = VELZ_VAR
     VELT1_VAR= VELX_VAR
     VELT2_VAR= VELY_VAR     
  end select

  ! Set the required array length
  sweepBegin = NGUARD-1
  sweepEnd   = n-NGUARD+1

  call hy_rhd_enthalpy (Vp, hp, sweepBegin,sweepEnd, n)  
  call hy_rhd_enthalpy (Vm, hm, sweepBegin,sweepEnd, n)

  ! Initialize array
  Vs = 0.

  ! Reset the required array length
  sweepBegin = NGUARD+1

  nfail = 0 
  do i = sweepBegin,sweepEnd

   !! solve Riemann problem between  qL(i-1/2) and qR(i-1/2)

    qR(:) = Vm(:,i)
    qL(:) = Vp(:,i-1)

    hR = hm(i)
    hL = hp(i-1)

    gR = 1.d0/sqrt(1.d0 - qR(VELX_VAR)*qR(VELX_VAR) &
                        - qR(VELY_VAR)*qR(VELY_VAR) &
                        - qR(VELZ_VAR)*qR(VELZ_VAR))

    gL = 1.d0/sqrt(1.d0 - qL(VELX_VAR)*qL(VELX_VAR) &
                        - qL(VELY_VAR)*qL(VELY_VAR) &
                        - qL(VELZ_VAR)*qL(VELZ_VAR))

    tauR = 1.d0/qR(DENS_VAR)
    tauL = 1.d0/qL(DENS_VAR)

    uxR = qR(VELN_VAR)
    uxL = qL(VELN_VAR)

    pR = qR(PRES_VAR)
    pL = qL(PRES_VAR)

    zeta_R = pR*tauR
    zeta_L = pL*tauL

    VR = tauR/gR
    VL = tauL/gL

    ! -----------------------------------------
    !   First guess is done outside the loop
    ! -----------------------------------------
    j2R = tauR*(hR*(gamma_r - 2.d0) + 1.d0)/(gamma_r*pR)
    j2L = tauL*(hL*(gamma_r - 2.d0) + 1.d0)/(gamma_r*pL)

    gR1    = sqrt(VR*VR + (1.d0 - uxR*uxR)*j2R)  !  ** RIGHT  
    zeta_R = VR*uxR + gR1

    gL1    = -sqrt(VL*VL + (1.d0 - uxL*uxL)*j2L) !  ** LEFT  
    zeta_L = VL*uxL + gL1

    duR = gR1/(hR*gR)
    duL = gL1/(hL*gL)

    p1 = uxL - uxR + duR*pR - duL*pL
    p1 = p1/(duR - duL)

    if (p1 .lt. 0.d0) then 
      p1 = 0.5d0*min(pL, pR)
    end if

    do iter = 1, MAX_ITER    ! ====== begin iteration loop =======
      call hy_rhd_shock (tauL, uxL, pL, gL, VL, hL, p1, uxL1, duL, zeta_L, -1)
      call hy_rhd_shock (tauR, uxR, pR, gR, VR, hR, p1, uxR1, duR, zeta_R,  1)

      dp = (uxR1 - uxL1)/(duL - duR)
      p1 = p1 + dp 
      if (abs(dp) .lt. 1.d-6*p1) exit
      if (p1 < 0.d0) then 
        p1 = 0.5d0*(p1 - dp)
      end if
    end do                   ! ====== end iteration loop ======

    if ((iter .ge. MAX_ITER) .or. (p1 .lt. 0.0)) then
      nfail           = nfail + 1
      failzone(nfail) = i
      Vs(:,i) = Vc(:,i)
      cycle
    end if
    u1 = 0.5d0*(uxR1 + uxL1)

    if (abs(u1) .lt. 1.d-10) then
      u1 = 0.d0
    end if

    Vstar(VELN_VAR) = u1
    Vstar(PRES_VAR) = p1

    ! ---------------------------------------------
    !      Sample solution on x/t = 0 axis
    ! ---------------------------------------------

    if (u1 .ge. 0.d0) then     !!   **** Left going Shock **** 

      gL1              = gL*hL/(gL*hL + (p1 - pL)*(VL + zeta_L*uxL))
      Vstar(VELT1_VAR) = qL(VELT1_VAR)*gL1
      Vstar(VELT2_VAR) = qL(VELT2_VAR)*gL1

      VL1        = VL - (u1 - uxL)*zeta_L
      gL1        = 1.d0/sqrt(1.d0 - Vstar(VELX_VAR)*Vstar(VELX_VAR) &
                                  - Vstar(VELY_VAR)*Vstar(VELY_VAR) &
                                  - Vstar(VELZ_VAR)*Vstar(VELZ_VAR))
          
      Vstar(DENS_VAR) = max(1.e-12, 1.0/(VL1*gL1))
      a          = VL/zeta_L + uxL       !   ** shock speed

      if (a .gt. 0.e0) then              !   ** region L
        Vs(:,i) = qL(:)
      else                               !   ** region L1
        Vs(:,i) = Vstar(:)
      end if
          
    else                     !!   **** Right going Shock **** 

      gR1        = gR*hR/(gR*hR + (p1 - pR)*(VR + zeta_R*uxR))
      Vstar(VELT1_VAR) = qR(VELT1_VAR)*gR1
      Vstar(VELT2_VAR) = qR(VELT2_VAR)*gR1

      gR1        = 1.e0/sqrt(1.e0 - Vstar(VELX_VAR)*Vstar(VELX_VAR) &
                                  - Vstar(VELY_VAR)*Vstar(VELY_VAR) &
                                  - Vstar(VELZ_VAR)*Vstar(VELZ_VAR))
      VR1        = VR - (u1 - uxR)*zeta_R
      Vstar(DENS_VAR) = max(1.e-12, 1.0/(VR1*gR1))
      a          = VR/zeta_R + uxR        !   ** shock speed

      if (a .gt. 0.e0) then               !   ** region R1
        Vs(:,i) = Vstar(:)
      else                                !   ** region R
        Vs(:,i) = qR(:)
      end if
    end if
  end do    ! ** close grid loop

! ---------------------------------------------------
!   Compute enthalpy, get fluxes, get maximum
!   characteristic speed
! ---------------------------------------------------
  call hy_rhd_flux(Vs, flux, Us, sweepBegin, sweepEnd, n, dir)
  call hy_rhd_enthalpy(Vs, hm, sweepBegin, sweepEnd, n)    
  call hy_rhd_maxChSpeed(Vs, hm, speed, sweepBegin, sweepEnd, n, dir)

! ---------------------------------------------------
!     Add artificial (Lapidus) Viscosity
! ---------------------------------------------------

  Uc = Vc
  call hy_rhd_primitiveToconserve(Uc, sweepBegin-1, sweepEnd, n)
 
  do i = sweepBegin, sweepEnd
    a = 0.1d0*max(0.e0, Vc(VELN_VAR, i-1) - Vc(VELN_VAR,i))
    flux(DENS_FLUX,i) = flux(DENS_FLUX,i) - a*(Uc(DENS_VAR,i) - Uc(DENS_VAR,i-1))
    flux(XMOM_FLUX,i) = flux(XMOM_FLUX,i) - a*(Uc(VELX_VAR,i) - Uc(VELX_VAR,i-1))
    flux(YMOM_FLUX,i) = flux(YMOM_FLUX,i) - a*(Uc(VELY_VAR,i) - Uc(VELY_VAR,i-1))
    flux(ZMOM_FLUX,i) = flux(ZMOM_FLUX,i) - a*(Uc(VELZ_VAR,i) - Uc(VELZ_VAR,i-1))
    flux(ENER_FLUX,i) = flux(ENER_FLUX,i) - a*(Uc(ENER_VAR,i) - Uc(ENER_VAR,i-1))
  end do

  ! --------------------------------------
  !  energy-pressure consistency relation
  ! --------------------------------------

  do i = sweepBegin, sweepEnd-1

    qR(DENS_VAR) = Uc(DENS_VAR,i) + dt/dx(i)*(flux(DENS_FLUX,i) - flux(DENS_FLUX,i+1))
    qR(ENER_VAR) = Uc(ENER_VAR,i) + dt/dx(i)*(flux(ENER_FLUX,i) - flux(ENER_FLUX,i+1))

    if (qR(DENS_VAR) .le. 0.0 .or. qR(ENER_VAR) .le. 0.0) then
      nfail = nfail + 1
      failzone(nfail) = i
    end if

!    qR(VELX_VAR) = Uc(VELX_VAR,i) + dt/dx(i)*(flux(VELX_VAR,i) - flux(VELX_VAR,i+1))
!    qR(VELY_VAR) = Uc(VELY_VAR,i) + dt/dx(i)*(flux(VELY_VAR,i) - flux(VELY_VAR,i+1))
!    qR(VELZ_VAR) = Uc(VELZ_VAR,i) + dt/dx(i)*(flux(VELZ_VAR,i) - flux(VELZ_VAR,i+1))
!    a = sqrt(qR(VELX_VAR)*qR(VELX_VAR) + qR(VELY_VAR)*qR(VELY_VAR) + qR(DENS_VAR)*qR(DENS_VAR))/qR(ENER_VAR)

!    if (a .ge. 1.e0) then
!      nfail = nfail + 1
!      failzone(nfail) = i
!    end if
  end do
     
  ! ----------------------------------------
  !   Replace failed Riemann problems with 
  !   a HLLE numerical flux
  ! ----------------------------------------

  if (nfail .gt. 0) then
    do iter = 1, nfail
      i = failzone(iter)
      !print *," ! Riemann solver has failed, fixing ...",i
      call hy_rhd_hlle(Vc, Vc, Vc, Flux, x, dx, dt, speed, vint, &
                       i - 1, i + 1, n, dir)
    end do
  end if

  return
end subroutine hy_rhd_riemann

