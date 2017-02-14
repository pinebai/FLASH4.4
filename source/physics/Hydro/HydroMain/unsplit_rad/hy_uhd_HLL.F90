!!****if* source/physics/Hydro/HydroMain/unsplit_rad/hy_uhd_HLL
!!
!! NAME
!!
!!  hy_uhd_HLL
!!
!! SYNOPSIS
!!
!!  call hy_uhd_HLL(integer(IN) :: dir,
!!             real(IN)    :: Vm(HY_VARINUMMAX),
!!             real(IN)    :: Vp(HY_VARINUMMAX),
!!             real(OUT)   :: Fstar(HY_VARINUM1),
!!             real(OUT)   :: speed,
!!             integer(OUT):: ierr)
!!
!! ARGUMENTS
!!
!!  dir    - a spatial direction for which the flux is being considered and computed
!!  Vm     - primitive variable for left state
!!            (DENS,VELX,VELY,VELZ,PRES,MAGX,MAGY,MAGZ + GAMC,GAME,EINT,TEMP)
!!  Vp     - primitive variable for right state
!!            (DENS,VELX,VELY,VELZ,PRES,MAGX,MAGY,MAGZ + GAMC,GAME,EINT,TEMP)
!!  Fstar  - computed flux data
!!            (includes face pressure at the end)
!!  speed  - fastest signal velocity to compute dt
!!  ierr   - a flag to check unphysical negative state (0 is ok; 1 is bad)
!!
!! DESCRIPTION
!! 
!!   This routine computes high-order Godunov fluxes based on the left and right Riemann states.
!!
!!   The HLL Riemann fan:
!!
!!            SL                  SR
!!             \                 /
!!              \               /
!!               \      U*     /
!!                \           /
!!                 \         /
!!                  \       /
!!           UL      \     /       UR
!!                    \   /
!!                     \ /
!!   --------------------------------------
!!
!! REFERENCES
!!
!!  * Harten, Lax and van Leer, SIAM  Review, 25(1):35--61, 1983
!!  * Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics, Springer, 1997
!!
!!***

Subroutine hy_uhd_HLL(dir,Vm,Vp,Fstar,speed,ierr)
  
  use hy_uhd_interface, ONLY : hy_uhd_prim2con,hy_uhd_prim2flx

  implicit none

#include "Flash.h"
#include "UHD.h"

  !! Arguments type declaration -----------
  integer, intent(IN) :: dir
  real, dimension(HY_VARINUMMAX), intent(IN)  :: Vm, Vp
  real, dimension(HY_VARINUM1),   intent(OUT) :: Fstar
  real,    intent(OUT) :: speed
  integer, intent(OUT) :: ierr
  !! --------------------------------------

  real :: SL,SR,cfL,cfR,aL2,aR2,velNL,velNR
  real :: magBL2,magBR2,magNL,magNR
  real, dimension(HY_VARINUM1) :: UL,UR,FL,FR
  real, parameter :: tiny=1.e-32
  integer :: hyEndVar,hyEndFlux


  ! Set index range depending on hydro or MHD
  ! default for hydro
  hyEndVar  = HY_ENER
  hyEndFlux = HY_P_FLUX
#ifdef FLASH_USM_MHD /* for USM-MHD */
  hyEndVar  = HY_MAGZ
#elif defined(FLASH_UGLM_MHD)
  hyEndVar  = HY_GLMP
  hyEndFlux = HY_GLMP_FLUX

  hyEndVar  = HY_MAGZ
  hyEndFlux = HY_MAGZ_FLUX
#endif



  ! Set no error to begin with
  ierr = 0

  ! Normal velocity
  velNL = Vm(HY_VELX+dir-1)
  velNR = Vp(HY_VELX+dir-1)

  ! Set sound speed
  aL2   = Vm(HY_GAMC)*Vm(HY_PRES)/Vm(HY_DENS)
  aR2   = Vp(HY_GAMC)*Vp(HY_PRES)/Vp(HY_DENS)

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD) /* compute additional MHD waves */
  magNL = Vm(HY_MAGX+dir-1)
  magNR = Vp(HY_MAGX+dir-1)
  magBL2= dot_product(Vm(HY_MAGX:HY_MAGZ),Vm(HY_MAGX:HY_MAGZ))/Vm(HY_DENS)
  magBR2= dot_product(Vp(HY_MAGX:HY_MAGZ),Vp(HY_MAGX:HY_MAGZ))/Vp(HY_DENS)
#endif

  ! Check unphysical negativity
  if ((Vm(HY_DENS) < tiny .and. Vm(HY_DENS) > 0.) .or. &
      (Vp(HY_DENS) < tiny .and. Vp(HY_DENS) > 0.) .or. &
      (Vm(HY_PRES) < tiny .and. Vm(HY_PRES) > 0.) .or. &
      (Vp(HY_PRES) < tiny .and. Vp(HY_PRES) > 0.)) then
     ! This could be vacuum limit. We return with zero flux.
     Fstar = 0.
     return
  elseif (aL2 < 0. .or. aR2 < 0.) then
     ierr = 1
     return
  else

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD) /* for MHD */
     cfL = sqrt(0.5*(aL2 + magBL2 + sqrt((aL2 + magBL2 )**2 - 4.*aL2*magNL*magNL/Vm(HY_DENS))))
     cfR = sqrt(0.5*(aR2 + magBR2 + sqrt((aR2 + magBR2 )**2 - 4.*aR2*magNR*magNR/Vp(HY_DENS))))
#else
     cfL = sqrt(aL2)
     cfR = sqrt(aR2)
#endif
  endif

  ! Get left/right going fastest wave speeds SL & SR for the left and right states
  ! by S. F. Davis, SIAM J. Sci. Stat, Comput., 9(1988) 445.
  ! Also see Miyoshi, Kusano, JCP, 208 (2005)
  SL = min(velNL - cfL, velNR - cfR)
  SR = max(velNL + cfL, velNR + cfR)

  ! Output maximum local wave speed for dt calculation
  !speed = abs(velNL)+0.5*(cfL+cfR)
  speed = max(abs(SL),abs(SR))

  ! Convert primitive variables to conservative variables
  call hy_uhd_prim2con(Vm(HY_DENS:HY_GAME),UL(HY_DENS:hyEndVar))
  call hy_uhd_prim2con(Vp(HY_DENS:HY_GAME),UR(HY_DENS:hyEndVar))
  UL(hyEndVar+1) = 0.0; UR(hyEndVar+1) = 0.0
  call hy_uhd_prim2flx(dir,Vm(HY_DENS:HY_VARINUM4),FL(F01DENS_FLUX:hyEndFlux))
  call hy_uhd_prim2flx(dir,Vp(HY_DENS:HY_VARINUM4),FR(F01DENS_FLUX:hyEndFlux))


  if (SL > 0.) then
     !Ustar(HY_DENS:HY_PRES) = UL(HY_DENS:HY_PRES)
     Fstar(F01DENS_FLUX:hyEndFlux) = FL(F01DENS_FLUX:hyEndFlux)
  elseif (SR < 0.) then
     !Ustar(HY_DENS:HY_PRES) = UR(HY_DENS:HY_PRES)
     Fstar(F01DENS_FLUX:hyEndFlux) = FR(F01DENS_FLUX:hyEndFlux)
  else !if ((SL <= 0.) .and. (SR >= 0.)) then
     !Ustar(HY_DENS:HY_PRES) = (SR*UR - SL*UL - FR + FL)/(SR - SL)
     Fstar(F01DENS_FLUX:hyEndFlux) = (  SR*FL(F01DENS_FLUX:hyEndFlux) &
                                      - SL*FR(F01DENS_FLUX:hyEndFlux) &
                                      + SR*SL*(UR(HY_DENS:hyEndVar+1)   &
                                             - UL(HY_DENS:hyEndVar+1))  &
                                      )/(SR - SL)
  endif

#ifdef FLASH_USM_MHD
  ! Enforce zero for corresponding magnetic field components
  Fstar(F06MAGX_FLUX+dir-1) = 0.
#endif

End Subroutine hy_uhd_HLL
