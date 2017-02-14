!!****if* source/physics/Hydro/HydroMain/unsplit/hy_uhd_Marquina
!!
!! NAME
!!
!!  hy_uhd_Marquina
!!
!! SYNOPSIS
!!
!!  call hy_uhd_Marquina( integer(IN) :: dir,
!!                   real(IN)    :: Vm(HY_VARINUMMAX),
!!                   real(IN)    :: Vp(HY_VARINUMMAX),
!!                   real(OUT)   :: Fstar(HY_VARINUM1),
!!                   real(OUT)   :: speed,
!!                   integer(OUT):: ierr)
!!
!! ARGUMENTS
!!
!!  dir    - a spatial direction for which the flux is being considered and computed
!!  Vm     - primitive variable for left state
!!            (DENS,VELX,VELY,VELZ,PRES,MAGX,MAGY,MAGZ + GAMC,GAME,EINT,TEMP)
!!  Vp     - primitive variable for right state
!!            (DENS,VELX,VELY,VELZ,PRES,MAGX,MAGY,MAGZ + GAMC,GAME,EINT,TEMP)
!!  Fstar  - computed flux data
!!            (includes face pressure at the end !!DEV: currently set to 0.0 !)
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

Subroutine hy_uhd_Marquina(dir,Vm,Vp,Fstar,speed,ierr)

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
  real, dimension(HY_VARINUMMAX), intent(IN) :: Vm, Vp
  real, dimension(HY_VARINUM1),  intent(OUT) :: Fstar
  real, intent(OUT) :: speed
  integer, intent(OUT) :: ierr
  !! --------------------------------------

  real,dimension(HY_VARINUM2) :: Vavg !DENS,VELX,VELY,VELZ,MAGX,MAGY,MAGZ,PRES + GAMC,GAME
  real,dimension(HY_VARINUM)   :: Um,Up,FLtemp,FRtemp,vecL,vecR
  real,dimension(HY_VARINUM1)  :: FL,FR
  real,dimension(HY_WAVENUM)   :: lambdaL,lambdaR
  real,dimension(HY_VARINUM,HY_WAVENUM):: reig,reigL,reigR,leig,leigL,leigR
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


  ! Godunov flux
  cons=.true.

#ifndef FLASH_USM_MHD
  ! Left state
  call hy_uhd_eigenParameters(Vm(HY_DENS:HY_GAME),dir,uN,cf)
  call hy_uhd_eigenValue(lambdaL,uN,cf)
  call hy_uhd_eigenVector(leigL,reigL,Vm(HY_DENS:HY_GAME),dir,cons,cf)

  ! Right state
  call hy_uhd_eigenParameters(Vp(HY_DENS:HY_GAME),dir,uN,cf)
  call hy_uhd_eigenValue(lambdaR,uN,cf)
  call hy_uhd_eigenVector(leigR,reigR,Vp(HY_DENS:HY_GAME),dir,cons,cf)
#else
  ! Avg state for MHD
  call hy_uhd_avgState(dir,Vm(HY_DENS:HY_EINT),Vp(HY_DENS:HY_EINT),Vavg(HY_DENS:HY_GAME))
  call hy_uhd_eigenParameters(Vavg(HY_DENS:HY_GAME),dir,uN,cf,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta)
  call hy_uhd_eigenVector(leigL,reigL,Vavg(HY_DENS:HY_GAME),dir,cons,cf,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta)

  leigR = leigL
  reigR = reigL

  ! Left state
  call hy_uhd_eigenParameters(Vm(HY_DENS:HY_GAME),dir,uN,cf,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta)
  call hy_uhd_eigenValue(lambdaL,uN,cf,C_alfn=ca,C_slow=cs)

  ! Right state
  call hy_uhd_eigenParameters(Vp(HY_DENS:HY_GAME),dir,uN,cf,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta)
  call hy_uhd_eigenValue(lambdaR,uN,cf,C_alfn=ca,C_slow=cs)
#endif


  ! Fluxes
  call hy_uhd_prim2con(Vm(HY_DENS:HY_GAME),Um(HY_DENS:hyEndVar))
  call hy_uhd_prim2con(Vp(HY_DENS:HY_GAME),Up(HY_DENS:hyEndVar))
  call hy_uhd_prim2flx(dir,Vm,FL(F01DENS_FLUX:hyEndFlux))
  call hy_uhd_prim2flx(dir,Vp,FR(F01DENS_FLUX:hyEndFlux))


 ! Output maximum local wave speed for dt calculation
  speed = max(abs(lambdaL(HY_FASTRGHT)),abs(lambdaR(HY_FASTRGHT)))

  FLtemp = 0.
  FRtemp = 0.

  do k=1,HY_WAVENUM
     if (lambdaL(k)*lambdaR(k) > 0.) then
        if (lambdaL(k) > 0.) then
           vecL(k)=dot_product(leigL(k,HY_DENS:hyEndVar),FL(F01DENS_FLUX:hyEndFlux-1)) &
                    + leigL(k,HY_XMOM+dir-1)           * FL(HY_P_FLUX)
           vecR(k)=0.
        else
           vecL(k)=0.
           vecR(k)=dot_product(leigR(k,HY_DENS:hyEndVar),FR(F01DENS_FLUX:hyEndFlux-1)) &
                    +          leigR(k,HY_XMOM+dir-1)  * FR(HY_P_FLUX)
        endif
     else
        vecL(k)=0.5*( dot_product(leigL(k,HY_DENS:hyEndVar),FL(F01DENS_FLUX:hyEndFlux-1))&
                                   +max(maxval(abs(lambdaL)),maxval(abs(lambdaR)))&
                                   *dot_product(leigL(k,HY_DENS:hyEndVar),Um(HY_DENS:hyEndVar)) &
                       +          leigL(k,HY_XMOM+dir-1)  * FL(HY_P_FLUX) )

        vecR(k)=0.5*( dot_product(leigR(k,HY_DENS:hyEndVar),FR(F01DENS_FLUX:hyEndFlux-1))&
                                   -max(maxval(abs(lambdaL)),maxval(abs(lambdaR)))&
                                   *dot_product(leigR(k,HY_DENS:hyEndVar),Up(HY_DENS:hyEndVar)) &
                       +          leigR(k,HY_XMOM+dir-1)  * FR(HY_P_FLUX) )
     endif
  enddo

  Fstar(HY_P_FLUX) = 0.0        !!DEV: What to do here instead?  - KW

  ! Left flux
  FLtemp(F01DENS_FLUX:hyEndFlux-1) = &
        vecL(HY_FASTLEFT)*reigL(HY_DENS:hyEndVar,HY_FASTLEFT)&
#ifdef FLASH_USM_MHD
       +vecL(HY_ALFNLEFT)*reigL(HY_DENS:hyEndVar,HY_ALFNLEFT)&
#endif
       +vecL(HY_SLOWLEFT)*reigL(HY_DENS:hyEndVar,HY_SLOWLEFT)&
       +vecL(HY_ENTROPY )*reigL(HY_DENS:hyEndVar,HY_ENTROPY) &
       +vecL(HY_SLOWRGHT)*reigL(HY_DENS:hyEndVar,HY_SLOWRGHT)&
#ifdef FLASH_USM_MHD
       +vecL(HY_ALFNRGHT)*reigL(HY_DENS:hyEndVar,HY_ALFNRGHT)&
#endif
       +vecL(HY_FASTRGHT)*reigL(HY_DENS:hyEndVar,HY_FASTRGHT)


  ! Right flux
  FRtemp(F01DENS_FLUX:hyEndFlux-1) = &
        vecR(HY_FASTLEFT)*reigR(HY_DENS:HY_ENER,HY_FASTLEFT)&
#ifdef FLASH_USM_MHD
       +vecR(HY_ALFNLEFT)*reigR(HY_DENS:hyEndVar,HY_ALFNLEFT)&
#endif
       +vecR(HY_SLOWLEFT)*reigR(HY_DENS:hyEndVar,HY_SLOWLEFT)&
       +vecR(HY_ENTROPY )*reigR(HY_DENS:hyEndVar,HY_ENTROPY) &
       +vecR(HY_SLOWRGHT)*reigR(HY_DENS:hyEndVar,HY_SLOWRGHT)&
#ifdef FLASH_USM_MHD
       +vecR(HY_ALFNRGHT)*reigR(HY_DENS:hyEndVar,HY_ALFNRGHT)&
#endif
       +vecR(HY_FASTRGHT)*reigR(HY_DENS:hyEndVar,HY_FASTRGHT)

  ! final flux in the star region
  Fstar(F01DENS_FLUX:hyEndFlux-1) = FLtemp(F01DENS_FLUX:hyEndFlux-1)&
                                   +FRtemp(F01DENS_FLUX:hyEndFlux-1)


#ifdef FLASH_USM_MHD
  ! Enforce zero for corresponding magnetic field components
  Fstar(F06MAGX_FLUX+dir-1) = 0.
#endif

End Subroutine hy_uhd_Marquina
