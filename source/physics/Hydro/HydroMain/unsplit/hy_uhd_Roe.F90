!!****if* source/physics/Hydro/HydroMain/unsplit/hy_uhd_Roe
!!
!! NAME
!!
!!  hy_uhd_Roe
!!
!! SYNOPSIS
!!
!!  hy_uhd_Roe( integer(IN) :: dir,
!!              real(IN)    :: Vm(HY_VARINUMMAX),
!!              real(IN)    :: Vp(HY_VARINUMMAX),
!!              real(OUT)   :: Fstar(HY_VARINUM1),
!!              real(OUT)   :: speed,
!!              integer(OUT):: ierr)
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
!! REFERENCES
!!
!!  * Roe, JCP, 43:357-372, 1981
!!
!!***

Subroutine hy_uhd_Roe(dir,Vm,Vp,Fstar,speed,ierr)

  use hy_uhd_interface, ONLY : hy_uhd_avgState,       &
                               hy_uhd_eigenParameters,&
                               hy_uhd_eigenValue,     &
                               hy_uhd_eigenVector,    &
                               hy_uhd_prim2con,       &
                               hy_uhd_prim2flx

  use Hydro_data,       ONLY : hy_entropy

  implicit none

#include "constants.h"
#include "Flash.h"
#include "UHD.h"

  !! Arguments type declaration -----------
  integer, intent(IN) :: dir
  real, dimension(HY_VARINUMMAX), intent(IN)  :: Vm, Vp
  real, dimension(HY_VARINUM1),   intent(OUT) :: Fstar
  real,    intent(OUT) :: speed
  integer, intent(OUT) :: ierr
  !! --------------------------------------

  real, dimension(HY_VARINUM2) :: Vavg !DENS,VELX,VELY,VELZ,PRES + GAMC,GAME
  real, dimension(HY_VARINUM)  :: Um,Up,vec
  real, dimension(HY_VARINUM1) :: sigF,FL,FR
  real, dimension(HY_WAVENUM)  :: lambda,lambdaL,lambdaR
  real, dimension(HY_VARINUM,HY_WAVENUM) :: reig,leig
  integer :: k
  logical :: cons
  real    :: uN,aL2,aR2,cf,cs,ca,as,af,bbN
  real, dimension(MDIM) :: beta
  real, parameter :: tiny=1.e-32
  integer :: hyEndVar,hyEndFlux


  ! Set index range depending on hydro or MHD
#ifndef FLASH_USM_MHD /* for hydro */
  hyEndVar  = HY_ENER
!  hyEndFlux = F05ENER_FLUX + 1
#else /* for MHD */
  hyEndVar  = HY_MAGZ
!  hyEndFlux = F08MAGZ_FLUX + 1
#endif
  hyEndFlux = HY_P_FLUX

  ierr = 0 ! Set no error

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
#ifdef FLASH_USM_MHD /* extra optional arguments for MHD */
       ,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta&
#endif
       )

  call hy_uhd_eigenValue&
       (lambda,uN,cf&
#ifdef FLASH_USM_MHD
       ,C_alfn=ca,C_slow=cs&
#endif
       )

  call hy_uhd_eigenVector&
       (leig,reig,Vavg(HY_DENS:HY_GAME),dir,cons,cf&
#ifdef FLASH_USM_MHD
       ,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta&
#endif
       )


  ! Left state
  call hy_uhd_prim2con(Vm(HY_DENS:HY_GAME),Um(HY_DENS:hyEndVar))
  call hy_uhd_prim2flx(dir,Vm,FL(F01DENS_FLUX:hyEndFlux))

  ! Right state
  call hy_uhd_prim2con(Vp(HY_DENS:HY_GAME),Up(HY_DENS:hyEndVar))
  call hy_uhd_prim2flx(dir,Vp,FR(F01DENS_FLUX:hyEndFlux))

  if (hy_entropy) then
     cons=.false.
     ! Entropy fix for low density region
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

     call hy_uhd_entropyFix(lambda,lambdaL,lambdaR)
  endif !end if (hy_entropy) then

  ! Output maximum local wave speed for dt calculation
  speed = abs(uN) + abs(cf)

  ! Initialize sigmasum sigF with zero
  sigF =0.
  do k=1,HY_WAVENUM
     vec(HY_DENS:hyEndVar)&
          = abs(lambda(k))*reig(HY_DENS:hyEndVar,k)&
           *dot_product(leig(HY_DENS:hyEndVar,k),Up(HY_DENS:hyEndVar)-Um(HY_DENS:hyEndVar))

     sigF(F01DENS_FLUX:hyEndFlux-1) = sigF(F01DENS_FLUX:hyEndFlux-1) + vec(HY_DENS:hyEndVar)
  enddo

  ! Godunov Flux
  Fstar(F01DENS_FLUX:hyEndFlux) = (  FR(F01DENS_FLUX:hyEndFlux) &
                                    +FL(F01DENS_FLUX:hyEndFlux) &
                                  -sigF(F01DENS_FLUX:hyEndFlux) &
                                  )*0.5

#ifdef FLASH_USM_MHD
  ! Enforce zero for corresponding magnetic field components
  Fstar(F06MAGX_FLUX+dir-1) = 0.
#endif

End Subroutine hy_uhd_Roe
