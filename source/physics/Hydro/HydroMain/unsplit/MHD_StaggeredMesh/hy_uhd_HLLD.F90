!!****if* source/physics/Hydro/HydroMain/unsplit/MHD_StaggeredMesh/hy_uhd_HLLD
!!
!! NAME
!!
!!  hy_uhd_HLLD
!!
!! SYNOPSIS
!!
!!  hy_uhd_HLLD( integer(IN) :: dir,
!!               real(IN)    :: Vm(HY_VARINUMMAX),
!!               real(IN)    :: Vp(HY_VARINUMMAX),
!!               real(OUT)   :: Fstar(HY_VARINUM1),
!!               real(OUT)   :: speed,
!!               integer(OUT):: ierr)
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
!!   The HLLD Riemann fan:
!!
!!            SL*       SM       SR*
!!   SL        \        |        /        SR      
!!     *        \       |       /        *
!!       *   UL* \ UL** | UR** / UR*   *
!!         *      \     |     /      *
!!           *     \    |    /     *
!!             *    \   |   /    *
!!           UL  *   \  |  /   *   UR
!!                 *  \ | /  *
!!                   * \|/ *
!!   --------------------------------------
!!
!! REFERENCE
!!
!!  * Miyoshi & Kusano, JCP, 208:315-344, 2005
!!
!!***

Subroutine hy_uhd_HLLD(dir,Vm,Vp,Fstar,speed,ierr)
  
  use hy_uhd_interface, ONLY : hy_uhd_prim2con,hy_uhd_prim2flx,hy_uhd_HLLC
  use Driver_interface, ONLY : Driver_abortFlash
  use hy_uhd_slopeLimiters, ONLY : signum
  use Hydro_data,       ONLY : hy_tiny

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

  real, parameter :: epsilon = 1.e-4
  real :: cfL,cfR,aL2,aR2,magBL2,magBR2
  real :: magnL,magnR,magtL,magtR,velnL,velnR
  real :: SM,SL,SR,SL2,SR2
  real :: dStarL,dStarR,prestL,prestR,pres,Bn_hll,denom
  real :: scrch1L,scrch1R,scrch2L,scrch2R
  real :: scrch3L,scrch3R,scrch4L,scrch4R,scrch5L,scrch5R
  real :: velxStarL,velyStarL,velzStarL,velxStarR,velyStarR,velzStarR
  real :: magxStarL,magyStarL,magzStarL,magxStarR,magyStarR,magzStarR
  real :: velxStar2,velyStar2,velzStar2,magxStar2,magyStar2,magzStar2
  real :: signumBn
  real, dimension(HY_VARINUM)  :: UL,UR,UCstarL,UCstarR,UCstar2L,UCstar2R
  real, dimension(HY_VARINUM1) :: FL,FR
  real, dimension(HY_MAGX:HY_MAGZ) :: Uhll(HY_MAGX:HY_MAGZ)
  logical :: degeneracyHLLD


  ! Set no error
  ierr = 0

  ! Set no degeneracy for MHD
  degeneracyHLLD=.false.

  ! Normal velocity
  velnL = Vm(HY_VELX+dir-1)
  velnR = Vp(HY_VELX+dir-1)

  ! Normal and transverse magnetic components
  magnL = Vm(HY_MAGX+dir-1)
  magnR = Vp(HY_MAGX+dir-1)

  ! Magnitude of B fields
  magBL2 = dot_product(Vm(HY_MAGX:HY_MAGZ),Vm(HY_MAGX:HY_MAGZ))
  magBR2 = dot_product(Vp(HY_MAGX:HY_MAGZ),Vp(HY_MAGX:HY_MAGZ))

  ! Total pressure
  prestL = Vm(HY_PRES) + 0.5*magBL2
  prestR = Vp(HY_PRES) + 0.5*magBR2

  ! Magnitude of transverse fields
  magtL = magBL2-magnL*magnL
  magtR = magBR2-magnR*magnR

  ! Normalize by density
  magBL2= magBL2/Vm(HY_DENS)
  magBR2= magBR2/Vp(HY_DENS)

  ! Set sound speed
  aL2 = Vm(HY_GAMC)*Vm(HY_PRES)/Vm(HY_DENS)
  aR2 = Vp(HY_GAMC)*Vp(HY_PRES)/Vp(HY_DENS)

  ! Fastest magnetoacoustic waves
  cfL = 0.5*(aL2+magBL2+sqrt((aL2+magBL2)**2-4.*aL2*magnL*magnL/Vm(HY_DENS)))
  cfR = 0.5*(aR2+magBR2+sqrt((aR2+magBR2)**2-4.*aR2*magnR*magnR/Vp(HY_DENS)))

  ! Check unphysical negativity
  if ((Vm(HY_DENS) < hy_tiny .and. Vm(HY_DENS) > 0.) .or. &
      (Vp(HY_DENS) < hy_tiny .and. Vp(HY_DENS) > 0.) .or. &
      (Vm(HY_PRES) < hy_tiny .and. Vm(HY_PRES) > 0.) .or. &
      (Vp(HY_PRES) < hy_tiny .and. Vp(HY_PRES) > 0.)) then
     ! This could be vacuum limit. We return with zero flux.
     Fstar = 0.
     return
  elseif (aL2 < 0. .or. aR2 < 0.) then
     ierr = 1
     return
  else
     cfL   = sqrt(cfL)
     cfR   = sqrt(cfR)
  endif


  ! Get left/right going fastest wave speeds SL & SR for the left and right states
  ! by S. F. Davis, SIAM J. Sci. Stat, Comput., 9(1988) 445.
  ! Also see Miyoshi, Kusano, JCP, 208 (2005)
  SL = min(velnL - cfL, velnR - cfR)
  SR = max(velnL + cfL, velnR + cfR)

  ! Output maximum local wave speed for dt calculation
  speed = max(abs(SL),abs(SR))

  ! Convert primitive variables to conservative variables
  call hy_uhd_prim2con(Vm(HY_DENS:HY_GAME),UL(HY_DENS:HY_MAGZ))
  call hy_uhd_prim2con(Vp(HY_DENS:HY_GAME),UR(HY_DENS:HY_MAGZ))
  call hy_uhd_prim2flx(dir,Vm(HY_DENS:HY_GAME),FL(F01DENS_FLUX:HY_P_FLUX))
  call hy_uhd_prim2flx(dir,Vp(HY_DENS:HY_GAME),FR(F01DENS_FLUX:HY_P_FLUX))


  ! Get magnetic field components of HLL states for Bn_hll
  if (SL > 0.) then
     Uhll(HY_MAGX:HY_MAGZ) = UL(HY_MAGX:HY_MAGZ)
  elseif ((SL <= 0.) .and. (SR >= 0.)) then
     Uhll(HY_MAGX:HY_MAGZ) = (  SR*UR(HY_MAGX:HY_MAGZ) &
                              - SL*UL(HY_MAGX:HY_MAGZ) &
                              )/(SR - SL)
  else
     Uhll(HY_MAGX:HY_MAGZ) = UR(HY_MAGX:HY_MAGZ)
  endif


  ! Calculate intermediate states ---------------------------------------------
  Bn_hll = Uhll(HY_PRES+dir) !=(SR*magnR-SL*magnL)/(SR-SL)


  !!***************************************
  !! (I)    UL* and UR* regions           *
  !!***************************************
  ! Normal velocity component and the middle wave SM
  ! SM = u*L = u**L = u*R = u**R
  SM =(Vp(HY_DENS)*velnR*(SR-velnR)-Vm(HY_DENS)*velnL*(SL-velnL)&
       -prestR+prestL-magnL*magnL+magnR*magnR)/&
      (Vp(HY_DENS)*(SR-velnR)-Vm(HY_DENS)*(SL-velnL))


  ! Convenient parameters
  scrch1L = SL-velnL
  scrch2L = SL-SM
  scrch3L = SM-velnL

  scrch1R = SR-velnR
  scrch2R = SR-SM
  scrch3R = SM-velnR


  ! Total pressure in the whole Riemann fan
  ! pres*L = pres*R = pres**L = pres**R = pres
  pres = scrch1R*Vp(HY_DENS)*prestL &
        -scrch1L*Vm(HY_DENS)*prestR &
        +Vm(HY_DENS)*Vp(HY_DENS)*scrch1R*scrch1L*(velnR-velnL)

  pres = pres/(scrch1R*Vp(HY_DENS)-scrch1L*Vm(HY_DENS))

  ! Densities in UL* and UR*
  dStarL = UL(HY_DENS)*scrch1L/scrch2L
  dStarR = UR(HY_DENS)*scrch1R/scrch2R

  SL2 = SM - abs(Bn_hll)/sqrt(dStarL) ! = SL*
  SR2 = SM + abs(Bn_hll)/sqrt(dStarR) ! = SR*

  ! Check if degeneracy happens: 
  ! This is the case of cf=ca, ca>a
  ! (1) transverse B=0 (note: in general, when trans B=0, cf=ca when ca>a; cs=a when ca<a), and
  ! (2) normal B .ne. 0, and strong to give ca>a.
  ! If this happens, we use HLLC (see PLUTO implementation as well)
  if ( (SL2-SL) < epsilon*(SM-SL) ) then
     degeneracyHLLD=.true.
  endif
  if ( (SR-SR2) < epsilon*(SR-SM) ) then
     degeneracyHLLD=.true.
  endif
  if (degeneracyHLLD) then
     call hy_uhd_HLLC(dir,Vm,Vp,Fstar,speed,ierr)
     return
  endif




  denom   = Vm(HY_DENS)*scrch1L*scrch2L-magnL*magnL
  scrch4L = scrch3L/denom
  scrch5L = (Vm(HY_DENS)*scrch1L*scrch1L-magnL*magnL)/denom

  denom   = Vp(HY_DENS)*scrch1R*scrch2R-magnR*magnR
  scrch4R = scrch3R/denom
  scrch5R = (Vp(HY_DENS)*scrch1R*scrch1R-magnR*magnR)/denom


  !! +++++++++++++++++++++++++++++++++++++++++++++!
  !!    Proceed calculating left star regions     !
  !! +++++++++++++++++++++++++++++++++++++++++++++!
  ! Left primitive variables
  select case (dir)
  case (DIR_X)
     magxStarL = Bn_hll
     velxStarL = SM
     magyStarL = Vm(HY_MAGY)*scrch5L
     magzStarL = Vm(HY_MAGZ)*scrch5L
     velyStarL = Vm(HY_VELY) - Vm(HY_MAGX)*Vm(HY_MAGY)*scrch4L
     velzStarL = Vm(HY_VELZ) - Vm(HY_MAGX)*Vm(HY_MAGZ)*scrch4L
  case (DIR_Y)
     magyStarL = Bn_hll
     velyStarL = SM
     magxStarL = Vm(HY_MAGX)*scrch5L
     magzStarL = Vm(HY_MAGZ)*scrch5L
     velxStarL = Vm(HY_VELX) - Vm(HY_MAGY)*Vm(HY_MAGX)*scrch4L
     velzStarL = Vm(HY_VELZ) - Vm(HY_MAGY)*Vm(HY_MAGZ)*scrch4L
  case (DIR_Z)
     magzStarL = Bn_hll
     velzStarL = SM
     magxStarL = Vm(HY_MAGX)*scrch5L
     magyStarL = Vm(HY_MAGY)*scrch5L
     velxStarL = Vm(HY_VELX) - Vm(HY_MAGZ)*Vm(HY_MAGX)*scrch4L
     velyStarL = Vm(HY_VELY) - Vm(HY_MAGZ)*Vm(HY_MAGY)*scrch4L
  end select

  ! Left conserved variables
  UCstarL(HY_DENS)= dStarL
  UCstarL(HY_VELX)= UCstarL(HY_DENS)*velxStarL
  UCstarL(HY_VELY)= UCstarL(HY_DENS)*velyStarL
  UCstarL(HY_VELZ)= UCstarL(HY_DENS)*velzStarL

  UCstarL(HY_MAGX)= magxStarL
  UCstarL(HY_MAGY)= magyStarL
  UCstarL(HY_MAGZ)= magzStarL
  UCstarL(HY_ENER)= scrch1L*UL(HY_ENER)-prestL*velnL+pres*SM+&
                    Bn_hll*(dot_product(Vm(HY_VELX:HY_VELZ),Vm(HY_MAGX:HY_MAGZ))&
                            -velxStarL*UCstarL(HY_MAGX)&
                            -velyStarL*UCstarL(HY_MAGY)&
                            -velzStarL*UCstarL(HY_MAGZ))
  UCstarL(HY_ENER) =  UCstarL(HY_ENER)/scrch2L



  !! +++++++++++++++++++++++++++++++++++++++++++++!
  !! (Id) Proceed calculating right star regions  !
  !! +++++++++++++++++++++++++++++++++++++++++++++!
  ! Right primitive variables
  select case (dir)
  case (DIR_X)
     magxStarR = Bn_hll
     velxStarR = SM
     magyStarR = Vp(HY_MAGY)*scrch5R
     magzStarR = Vp(HY_MAGZ)*scrch5R
     velyStarR = Vp(HY_VELY) - Vp(HY_MAGX)*Vp(HY_MAGY)*scrch4R
     velzStarR = Vp(HY_VELZ) - Vp(HY_MAGX)*Vp(HY_MAGZ)*scrch4R
  case (DIR_Y)
     magyStarR = Bn_hll
     velyStarR = SM
     magxStarR = Vp(HY_MAGX)*scrch5R
     magzStarR = Vp(HY_MAGZ)*scrch5R
     velxStarR = Vp(HY_VELX) - Vp(HY_MAGY)*Vp(HY_MAGX)*scrch4R
     velzStarR = Vp(HY_VELZ) - Vp(HY_MAGY)*Vp(HY_MAGZ)*scrch4R
  case (DIR_Z)
     magzStarR = Bn_hll
     velzStarR = SM
     magxStarR = Vp(HY_MAGX)*scrch5R
     magyStarR = Vp(HY_MAGY)*scrch5R
     velxStarR = Vp(HY_VELX) - Vp(HY_MAGZ)*Vp(HY_MAGX)*scrch4R
     velyStarR = Vp(HY_VELY) - Vp(HY_MAGZ)*Vp(HY_MAGY)*scrch4R
  end select

  ! Right conserved variables
  UCstarR(HY_DENS)= dStarR
  UCstarR(HY_VELX)= UCstarR(HY_DENS)*velxStarR
  UCstarR(HY_VELY)= UCstarR(HY_DENS)*velyStarR
  UCstarR(HY_VELZ)= UCstarR(HY_DENS)*velzStarR

  UCstarR(HY_MAGX)= magxStarR
  UCstarR(HY_MAGY)= magyStarR
  UCstarR(HY_MAGZ)= magzStarR
  UCstarR(HY_ENER)= scrch1R*UR(HY_ENER)-prestR*velnR+pres*SM+&
                    Bn_hll*(dot_product(Vp(HY_VELX:HY_VELZ),Vp(HY_MAGX:HY_MAGZ))&
                           -velxStarR*UCstarR(HY_MAGX)&
                           -velyStarR*UCstarR(HY_MAGY)&
                           -velzStarR*UCstarR(HY_MAGZ))
  UCstarR(HY_ENER) = UCstarR(HY_ENER)/scrch2R
  !! Done with calculating UL* and UR* regions !!



  !!***************************************
  !! (II)    UL** and UR** regions        *
  !!***************************************
  ! First calculate SL* and SR*
  SL2 = SM - abs(Bn_hll)/sqrt(UCstarL(HY_DENS)) ! = SL*
  SR2 = SM + abs(Bn_hll)/sqrt(UCstarR(HY_DENS)) ! = SR*

  ! Densities
  UCstar2L(HY_DENS) = UCstarL(HY_DENS)
  UCstar2R(HY_DENS) = UCstarR(HY_DENS)

  scrch1L = sqrt(UCstarL(HY_DENS))
  scrch1R = sqrt(UCstarR(HY_DENS))
  scrch2L = 1./(scrch1L + scrch1R)
  scrch2R = scrch2L

  signumBn = signum(Bn_hll)

  select case (dir)
  case (DIR_X)
     ! Left primitive variables
     velxStar2 = SM
     velyStar2 = (scrch1L*velyStarL+scrch1R*velyStarR&
                +(UCstarR(HY_MAGY)-UCstarL(HY_MAGY))*signumBn)*scrch2L
     velzStar2 = (scrch1L*velzStarL+scrch1R*velzStarR&
                +(UCstarR(HY_MAGZ)-UCstarL(HY_MAGZ))*signumBn)*scrch2L

     magxStar2 = Bn_hll
     magyStar2 = (scrch1L*magyStarR+scrch1R*magyStarL&
                 +scrch1L*scrch1R*(velyStarR-velyStarL)*signumBn)&
                 *scrch2L
     magzStar2 = (scrch1L*magzStarR+scrch1R*magzStarL&
                 +scrch1L*scrch1R*(velzStarR-velzStarL)*signumBn)&
                 *scrch2L

  case (DIR_Y)
     ! Left primitive variables
     velxStar2 = (scrch1L*velxStarL+scrch1R*velxStarR&
                +(UCstarR(HY_MAGX)-UCstarL(HY_MAGX))*signumBn)*scrch2L
     velyStar2 = SM
     velzStar2 = (scrch1L*velzStarL+scrch1R*velzStarR&
                +(UCstarR(HY_MAGZ)-UCstarL(HY_MAGZ))*signumBn)*scrch2L

     magxStar2 = (scrch1L*magxStarR+scrch1R*magxStarL&
                 +scrch1L*scrch1R*(velxStarR-velxStarL)*signumBn)&
                 *scrch2L

     magyStar2 = Bn_hll
     magzStar2 = (scrch1L*magzStarR+scrch1R*magzStarL&
                 +scrch1L*scrch1R*(velzStarR-velzStarL)*signumBn)&
                 *scrch2L

  case (DIR_Z)
     ! Left primitive variables
     velxStar2 = (scrch1L*velxStarL+scrch1R*velxStarR&
                +(UCstarR(HY_MAGX)-UCstarL(HY_MAGX))*signumBn)*scrch2L
     velyStar2 = (scrch1L*velyStarL+scrch1R*velyStarR&
                +(UCstarR(HY_MAGY)-UCstarL(HY_MAGY))*signumBn)*scrch2L
     velzStar2 = SM

     magxStar2 = (scrch1L*magxStarR+scrch1R*magxStarL&
                 +scrch1L*scrch1R*(velxStarR-velxStarL)*signumBn)&
                 *scrch2L
     magyStar2 = (scrch1L*magyStarR+scrch1R*magyStarL&
                 +scrch1L*scrch1R*(velyStarR-velyStarL)*signumBn)&
                 *scrch2L
     magzStar2 = Bn_hll

  end select

  ! Left conservative variables
  UCstar2L(HY_VELX) = UCstar2L(HY_DENS)*velxStar2
  UCstar2L(HY_VELY) = UCstar2L(HY_DENS)*velyStar2
  UCstar2L(HY_VELZ) = UCstar2L(HY_DENS)*velzStar2

  UCstar2L(HY_MAGX) = magxStar2
  UCstar2L(HY_MAGY) = magyStar2
  UCstar2L(HY_MAGZ) = magzStar2
  UCstar2L(HY_ENER) = UCstarL(HY_ENER)-sqrt(UCstarL(HY_DENS))*signumBn*&
                      (velxStarL*UCstarL (HY_MAGX)&
                      +velyStarL*UCstarL (HY_MAGY)&
                      +velzStarL*UCstarL (HY_MAGZ)&
                      -velxStar2*UCstar2L(HY_MAGX)&
                      -velyStar2*UCstar2L(HY_MAGY)&
                      -velzStar2*UCstar2L(HY_MAGZ))

  ! Right conservative variables
  UCstar2R(HY_VELX) = UCstar2R(HY_DENS)*velxStar2
  UCstar2R(HY_VELY) = UCstar2R(HY_DENS)*velyStar2
  UCstar2R(HY_VELZ) = UCstar2R(HY_DENS)*velzStar2

  UCstar2R(HY_MAGX) = magxStar2
  UCstar2R(HY_MAGY) = magyStar2
  UCstar2R(HY_MAGZ) = magzStar2
  UCstar2R(HY_ENER) = UCstarR(HY_ENER)+sqrt(UCstarR(HY_DENS))*signumBn*&
                      (velxStarR*UCstarR (HY_MAGX)&
                      +velyStarR*UCstarR (HY_MAGY)&
                      +velzStarR*UCstarR (HY_MAGZ)&
                      -velxStar2*UCstar2R(HY_MAGX)&
                      -velyStar2*UCstar2R(HY_MAGY)&
                      -velzStar2*UCstar2R(HY_MAGZ))
  !! END of calculating all HLLD states !!

  !!***************************************
  !! (III) HLLD fluxes                    *
  !!***************************************
  if (SL >= 0.) then
     Fstar(F01DENS_FLUX:HY_P_FLUX) = FL(F01DENS_FLUX:HY_P_FLUX)

  elseif ((SL < 0.) .and. (SL2 >= 0.)) then
     Fstar(F01DENS_FLUX:F08MAGZ_FLUX) = FL(F01DENS_FLUX:F08MAGZ_FLUX) &
          + SL*(UCstarL(HY_DENS:HY_MAGZ) - UL(HY_DENS:HY_MAGZ))
     Fstar(HY_P_FLUX) = FL(HY_P_FLUX)

  elseif ((SL2 < 0.) .and. (SM >= 0.)) then
     Fstar(F01DENS_FLUX:F08MAGZ_FLUX) = FL(F01DENS_FLUX:F08MAGZ_FLUX) &
          + SL2*(UCstar2L(HY_DENS:HY_MAGZ) - UCstarL(HY_DENS:HY_MAGZ))&
          +  SL*( UCstarL(HY_DENS:HY_MAGZ) - UL(HY_DENS:HY_MAGZ))
     Fstar(HY_P_FLUX) = FL(HY_P_FLUX)

  elseif ((SM < 0.) .and. (SR2 >= 0.)) then
     Fstar(F01DENS_FLUX:F08MAGZ_FLUX) = FR(F01DENS_FLUX:F08MAGZ_FLUX) &
          + SR2*(UCstar2R(HY_DENS:HY_MAGZ) - UCstarR(HY_DENS:HY_MAGZ))&
          +  SR*( UCstarR(HY_DENS:HY_MAGZ) - UR(HY_DENS:HY_MAGZ))
     Fstar(HY_P_FLUX) = FR(HY_P_FLUX)

  elseif ((SR2 < 0.) .and. (SR >= 0.)) then
     Fstar(F01DENS_FLUX:F08MAGZ_FLUX) = FR(F01DENS_FLUX:F08MAGZ_FLUX) &
          + SR*(UCstarR(HY_DENS:HY_MAGZ) - UR(HY_DENS:HY_MAGZ))
     Fstar(HY_P_FLUX) = FR(HY_P_FLUX)

  else
     Fstar(F01DENS_FLUX:HY_P_FLUX) = FR(F01DENS_FLUX:HY_P_FLUX)

  endif

  !!!***************************************
  !! (IV) Enforce zero for corresponding   *
  !!      magnetic field components        *
  !!!***************************************
  Fstar(F06MAGX_FLUX+dir-1) = 0.


End Subroutine hy_uhd_HLLD
