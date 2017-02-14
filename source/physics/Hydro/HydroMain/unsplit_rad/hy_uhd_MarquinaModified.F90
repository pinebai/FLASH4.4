!!****if* source/physics/Hydro/HydroMain/unsplit_rad/hy_uhd_MarquinaModified
!!
!! NAME
!!
!!  hy_uhd_MarquinaModified
!!
!! SYNOPSIS
!!
!!  call hy_uhd_MarquinaModified( integer(IN) :: dir,
!!                           real(IN)    :: Vm(HY_VARINUMMAX),
!!                           real(IN)    :: Vp(HY_VARINUMMAX),
!!                           real(OUT)   :: Fstar(HY_VARINUM1),
!!                           real(OUT)   :: speed,
!!                           integer(OUT):: ierr)
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
!!  * Donat and Marquina, JCP, 125:42-58, 1996
!!  * Stiriba and Donat, Computers & Mathematics with Applications, 46:719-739, 2003
!!  * Leveque, Mihalas, Dorfi, and Muller, Computational Methods for Astrophysical Flows, Springer, 1997
!!  * Cunningham et al., AstroBear
!!
!!***

Subroutine hy_uhd_MarquinaModified(dir,Vm,Vp,Fstar,speed,ierr)

  use hy_uhd_interface, ONLY : hy_uhd_eigenParameters,&
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
  real, dimension(HY_VARINUM1),    intent(OUT) :: Fstar
  real,    intent(OUT) :: speed
  integer, intent(OUT) :: ierr
  !! --------------------------------------

  real,dimension(HY_VARINUM)  :: Um,Up,vec
  real,dimension(HY_VARINUM1) :: sigF,FL,FR
  real,dimension(HY_WAVENUM)  :: lambda,lambdaL,lambdaR,lambL,lambR
  real,dimension(HY_VARINUM,HY_WAVENUM):: leig,leigL,leigR,reig,reigL,reigR
  integer :: k
  logical :: cons
  real    :: uN,aL2,aR2,cf,cs,ca,as,af,bbN
  real, dimension(MDIM) :: beta
  real, parameter :: tiny=1.e-32
  integer :: hyEndVar,hyEndFlux


  ! Set index range depending on hydro or MHD
#ifndef FLASH_USM_MHD /* for hydro */
  hyEndVar  = HY_ENER
  hyEndFlux = F05ENER_FLUX + 1
#else /* for MHD */
  hyEndVar  = HY_MAGZ
  hyEndFlux = F08MAGZ_FLUX + 1
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


  ! Second order Godunov flux
  cons=.true.

  ! Left state
  call hy_uhd_eigenParameters&
       (Vm(HY_DENS:HY_GAME),dir,uN,cf&
#ifdef FLASH_USM_MHD /* extra optional arguments for MHD */
       ,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta&
#endif
       )

  call hy_uhd_eigenValue&
       (lambdaL,uN,cf&
#ifdef FLASH_USM_MHD
       ,C_alfn=ca,C_slow=cs&
#endif
       )

  call hy_uhd_eigenVector&
       (leigL,reigL,Vm(HY_DENS:HY_GAME),dir,cons,cf&
#ifdef FLASH_USM_MHD
       ,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta&
#endif
       )

  call hy_uhd_prim2con(Vm(HY_DENS:HY_GAME),Um(HY_DENS:hyEndVar))
  call hy_uhd_prim2flx(dir,Vm(HY_DENS:HY_VARINUM4),FL(F01DENS_FLUX:hyEndFlux))

  ! Right state
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

  call hy_uhd_eigenVector&
       (leigR,reigR,Vp(HY_DENS:HY_GAME),dir,cons,cf&
#ifdef FLASH_USM_MHD
       ,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta&
#endif
       )

  call hy_uhd_prim2con(Vp(HY_DENS:HY_GAME),Up(HY_DENS:hyEndVar))
  call hy_uhd_prim2flx(dir,Vp(HY_DENS:HY_VARINUM4),FR(F01DENS_FLUX:hyEndFlux))


  ! Output maximum local wave speed for dt calculation
  speed = abs(uN) + abs(cf)

  sigF =0.

  do k=1,HY_WAVENUM
     vec(HY_DENS:hyEndVar)&
          = 0.5*(abs(lambdaL(k))+abs(lambdaR(k)))*& ! this give circular symmetries
          ( reigR(HY_DENS:hyEndVar,k)*dot_product(leigR(HY_DENS:hyEndVar,k),Up(HY_DENS:hyEndVar)) &
           -reigL(HY_DENS:hyEndVar,k)*dot_product(leigL(HY_DENS:hyEndVar,k),Um(HY_DENS:hyEndVar)) )

     sigF(F01DENS_FLUX:hyEndFlux-1) = sigF(F01DENS_FLUX:hyEndFlux-1) + vec(HY_DENS:hyEndVar)
  enddo

  ! Godunov Flux
  Fstar(F01DENS_FLUX:hyEndFlux) = (   FR(F01DENS_FLUX:hyEndFlux) &
                                     +FL(F01DENS_FLUX:hyEndFlux) &
                                   -sigF(F01DENS_FLUX:hyEndFlux) &
                                   )*0.5


#ifdef FLASH_USM_MHD
  ! Enforce zero for corresponding magnetic field components
  Fstar(F06MAGX_FLUX+dir-1) = 0.
#endif

End Subroutine hy_uhd_MarquinaModified
