!!****if* source/physics/Hydro/HydroMain/unsplit_rad/hy_uhd_LLF
!!
!! NAME
!!
!!  hy_uhd_LLF
!!
!! SYNOPSIS
!!
!!  hy_uhd_LLF(integer(IN) :: dir,
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
!!  speed  - fastest signal velocity to compute dt
!!  ierr   - a flag to check unphysical negative state (0 is ok; 1 is bad)
!!
!! DESCRIPTION
!! 
!!   This routine computes local Lax-Friedrichs fluxes that are conservative and consistent at intercells.
!!   This flux has the maximum amount of dissipation allowed by the stability condition and
!!   exhibits odd-even decoupling. (See e.g., Laney, Computational Gasdynamics, pp.315)
!!
!! REFERENCE
!!
!!  * Laney, Computational Gasdynamics, Cambridge University Press
!!
!!***

Subroutine hy_uhd_LLF(dir,Vm,Vp,Fstar,speed,ierr)

  use hy_uhd_interface, ONLY : hy_uhd_avgState,       &
                               hy_uhd_prim2con,       &
                               hy_uhd_prim2flx,       &
                               hy_uhd_eigenParameters,&
                               hy_uhd_eigenValue

  implicit none

#include "Flash.h"
#include "constants.h"
#include "UHD.h"

  !! Arguments type declaration -----------
  integer, intent(IN) :: dir
  real, dimension(HY_VARINUMMAX), intent(IN)  :: Vm, Vp
  real, dimension(HY_VARINUM1),   intent(OUT) :: Fstar
  real,    intent(OUT) :: speed
  integer, intent(OUT) :: ierr
  !! --------------------------------------

  real,dimension(HY_VARINUM2) :: Vavg !DENS,VELX,VELY,VELZ,PRES + GAMC,GAME
  real,dimension(HY_VARINUM1) :: Um,Up,FL,FR
  real,dimension(HY_WAVENUM)  :: lambda,lambdaL,lambdaR
  logical :: cons
  real    :: uN,aL2,aR2,cf,cs,ca,as,af,bbN
  real, dimension(MDIM) :: beta
  real, parameter :: tiny=1.e-32
  integer :: hyEndVar,hyEndFlux

  ! Set index range depending on hydro or MHD
#ifndef FLASH_USM_MHD /* for hydro */
  hyEndVar  = HY_ENER
  hyEndFlux = F05ENER_FLUX+1
#else /* for MHD */
  hyEndVar  = HY_MAGZ
  hyEndFlux = F08MAGZ_FLUX+1
#endif

  ! Set no error
  ierr = 0

  ! Set sound speed
  aL2 = Vm(HY_GAMC)*Vm(HY_PRES)/Vm(HY_DENS)
  aR2 = Vp(HY_GAMC)*Vp(HY_PRES)/Vp(HY_DENS)

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
  endif


  ! Godunov flux
  cons=.true.

  ! Averaged state
  call hy_uhd_avgState&
       (dir,Vm(HY_DENS:HY_EINT),Vp(HY_DENS:HY_EINT),Vavg(HY_DENS:HY_GAME))

  call hy_uhd_eigenParameters&
       (Vavg(HY_DENS:HY_GAME),dir,uN,cf&
#ifdef FLASH_USM_MHD
       ,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta&
#endif
       )

  call hy_uhd_eigenValue&
       (lambda,uN,cf&
#ifdef FLASH_USM_MHD
       ,C_alfn=ca,C_slow=cs&
#endif
       )

  ! Left state
  call hy_uhd_prim2con(Vm(HY_DENS:HY_GAME),Um(HY_DENS:hyEndVar))
  call hy_uhd_prim2flx(dir,Vm(HY_DENS:HY_VARINUM4),FL(F01DENS_FLUX:hyEndFlux))
  call hy_uhd_eigenParameters&
       (Vm(HY_DENS:HY_GAME),dir,uN,cf&
#ifdef FLASH_USM_MHD
       ,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta&
#endif
       )

  call hy_uhd_eigenValue&
       (lambdaL,uN,cf&
#ifdef FLASH_USM_MHD
       ,C_alfn=ca,C_slow=cs&
#endif
       )

  ! Right state
  call hy_uhd_prim2con(Vp(HY_DENS:HY_GAME),Up(HY_DENS:hyEndVar))
  call hy_uhd_prim2flx(dir,Vp(HY_DENS:HY_VARINUM4),FR(F01DENS_FLUX:hyEndFlux))
  call hy_uhd_eigenParameters&
       (Vp(HY_DENS:HY_GAME),dir,uN,cf&
#ifdef FLASH_USM_MHD
       ,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta&
#endif
       )
  call hy_uhd_eigenValue&
       (lambdaR,uN,cf&
#ifdef FLASH_USM_MHD
       ,C_alfn=ca,C_slow=cs&
#endif
       )


  ! Output maximum local wave speed for dt calculation
  ! speed = lambda(HY_FASTRGHT)
  speed = abs(uN) + abs(cf)

  Um(hyEndVar+1) = 0.0; Up(hyEndVar+1) = 0.0

  !! local Lax-Friedrichs Flux
  Fstar(F01DENS_FLUX:hyEndFlux)  = ( FR(F01DENS_FLUX:hyEndFlux) &
                                    +FL(F01DENS_FLUX:hyEndFlux) &
                                    -max(maxval(abs(lambda )),&
                                         maxval(abs(lambdaL)),&
                                         maxval(abs(lambdaR)) &
                                         )*(Up(HY_DENS:hyEndVar+1)-Um(HY_DENS:hyEndVar+1))&
                                   )*0.5


#ifdef FLASH_USM_MHD
  ! Enforce zero for corresponding magnetic field components
  Fstar(F06MAGX_FLUX+dir-1) = 0.
#endif

End Subroutine hy_uhd_LLF
