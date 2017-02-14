!!****if* source/physics/Hydro/HydroMain/unsplit_rad/hy_uhd_HLLC
!!
!! NAME
!!
!!  hy_uhd_HLLC
!!
!! SYNOPSIS
!!
!!  hy_uhd_HLLC( integer(IN) :: dir,
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
!!   The HLLC Riemann fan:
!!
!!            SL      SM=qStar    SR
!!             \        |        /
!!              \       |       /
!!               \  UL* | UR*  /
!!                \     |     /
!!                 \    |    /
!!                  \   |   /
!!           UL      \  |  /       UR
!!                    \ | /
!!                     \|/
!!   --------------------------------------
!!
!! REFERENCES
!!
!!  * S. Li, JCP, 203:344-357, 2005
!!  * Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics, Springer, 1997
!!
!!***

Subroutine hy_uhd_HLLC(dir,Vm,Vp,Fstar,speed,ierr)

  use hy_uhd_interface, ONLY : hy_uhd_prim2con,hy_uhd_prim2flx
  use Driver_interface, ONLY : Driver_abortFlash
  use Hydro_data,       ONLY : hy_tiny, fP => hy_fPresInMomFlux

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

  real :: SL,SR,cfL,cfR,aL2,aR2,velNL,velNR
  real :: dStarL,dStarR,totalPresL,totalPresR
  real :: BxStar,ByStar,BzStar,Bn_hll,pStar,qStar
  real :: denomL,denomR,numerL,numerR
  real, dimension(HY_VARINUM)  :: UL,UR,Uhll,UCstarR,UCstarL
  real, dimension(HY_VARINUM1) :: FL,FR
  real :: magBL2,magBR2,magNL,magNR
  integer :: hyEndVar,hyEndFlux

  ! Set index range depending on hydro or MHD
  hyEndVar  = HY_ENER
  hyEndFlux = HY_ENER_FLUX+1
#ifdef FLASH_USM_MHD /* for USM-MHD */
  hyEndVar  = HY_MAGZ
  hyEndFlux = HY_MAGZ_FLUX+1
#elif defined(FLASH_UGLM_MHD)
  hyEndVar  = HY_GLMP
  hyEndFlux = HY_GLMP_FLUX
#endif


  ! Set no error to begin with
  ierr = 0

  ! Normal velocity
  velNL = Vm(HY_VELX+dir-1)
  velNR = Vp(HY_VELX+dir-1)

  ! Set sound speed
  aL2   = Vm(HY_GAMC)*Vm(HY_PRES)/Vm(HY_DENS)
  aR2   = Vp(HY_GAMC)*Vp(HY_PRES)/Vp(HY_DENS)


  ! Set zero magnetic quantities by default for hydro
  magNL = 0.
  magNR = 0.

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD) /* compute additional MHD waves */
  magNL = Vm(HY_MAGX+dir-1)
  magNR = Vp(HY_MAGX+dir-1)
  magBL2= dot_product(Vm(HY_MAGX:HY_MAGZ),Vm(HY_MAGX:HY_MAGZ))/Vm(HY_DENS)
  magBR2= dot_product(Vp(HY_MAGX:HY_MAGZ),Vp(HY_MAGX:HY_MAGZ))/Vp(HY_DENS)
#endif

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
#ifdef FLASH_UHD_HYDRO /* for hydro */
     cfL = sqrt(aL2)
     cfR = sqrt(aR2)
#else /* for MHD */
     cfL = sqrt(0.5*(aL2 + magBL2 + sqrt((aL2 + magBL2 )**2 - 4.*aL2*magNL*magNL/Vm(HY_DENS))))
     cfR = sqrt(0.5*(aR2 + magBR2 + sqrt((aR2 + magBR2 )**2 - 4.*aR2*magNR*magNR/Vp(HY_DENS))))
#endif
  endif

  ! Get left/right going fastest wave speeds SL & SR for the left and right states
  ! by S. F. Davis, SIAM J. Sci. Stat, Comput., 9(1988) 445.
  ! Also see Miyoshi, Kusano, JCP, 208 (2005)
  SL = min(velNL - cfL, velNR - cfR)
  SR = max(velNL + cfL, velNR + cfR)

  ! Output maximum local wave speed for dt calculation
  speed = max(abs(SL),abs(SR))

  ! Total pressure
  totalPresL = Vm(HY_PRES)
  totalPresR = Vp(HY_PRES)
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD) /* for MHD */
  totalPresL = totalPresL + 0.5*dot_product(Vm(HY_MAGX:HY_MAGZ),Vm(HY_MAGX:HY_MAGZ))
  totalPresR = totalPresR + 0.5*dot_product(Vp(HY_MAGX:HY_MAGZ),Vp(HY_MAGX:HY_MAGZ))
#endif


  ! Convert primitive variables to conservative variables
  call hy_uhd_prim2con(Vm(HY_DENS:HY_GAME),UL(HY_DENS:hyEndVar))
  call hy_uhd_prim2con(Vp(HY_DENS:HY_GAME),UR(HY_DENS:hyEndVar))
  call hy_uhd_prim2flx(dir,Vm(HY_DENS:HY_VARINUM4),FL(F01DENS_FLUX:hyEndFlux))
  call hy_uhd_prim2flx(dir,Vp(HY_DENS:HY_VARINUM4),FR(F01DENS_FLUX:hyEndFlux))


  ! Get HLL states for later use
  if (SL > 0.) then
     Uhll(HY_DENS:hyEndVar) = UL(HY_DENS:hyEndVar)
  elseif ((SL <= 0.) .and. (SR >= 0.)) then
     Uhll(HY_DENS:hyEndVar) = &
          ( SR*UR(HY_DENS:hyEndVar) &
           -SL*UL(HY_DENS:hyEndVar) &
             - FR(F01DENS_FLUX:hyEndFlux-1) &
             + FL(F01DENS_FLUX:hyEndFlux-1)&
           )/(SR - SL)
     Uhll(HY_DENS+dir) = Uhll(HY_DENS+dir) + &
          (  totalPresL - totalPresR  )/(SR - SL) * (1.0-fP)
  else
     Uhll(HY_DENS:hyEndVar) = UR(HY_DENS:hyEndVar)
  endif

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD) /* for MHD */
  ! Calculate intermediate states --------------
  Bn_hll = Uhll(HY_PRES+dir) !=(SR*magNR-SL*magNL)/(SR-SL)
  BxStar = Uhll(HY_MAGX)     !BxStarL = BxStarR = BxHLL
  ByStar = Uhll(HY_MAGY)     !ByStarL = ByStarR = ByHLL
  BzStar = Uhll(HY_MAGZ)     !BzStarL = BzStarR = BzHLL
#endif

  ! (1) Normal velocity component
  ! qStarL = qStarR = qStar
  qStar =( Vp(HY_DENS)*velNR*(SR-velNR) &
          -Vm(HY_DENS)*velNL*(SL-velNL) &
          +totalPresL - totalPresR      &
          -magNL**2   + magNR**2 )
  qStar = qStar/(Vp(HY_DENS)*(SR-velNR) - Vm(HY_DENS)*(SL-velNL))

  ! Convenient parameters
  numerL = SL-velNL
  denomL = SL-qStar
  numerR = SR-velNR
  denomR = SR-qStar

  ! (2) Total pressure in the star region
  ! pStarL = pStarR = pStar
  pStar = Vm(HY_DENS)*numerL*(qStar-velNL)+totalPresL
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD) /* for MHD */
  pStar = pStar - magNL**2 + Bn_hll**2
#endif

  ! (3) Density
  dStarL = UL(HY_DENS)*numerL/denomL
  dStarR = UR(HY_DENS)*numerR/denomR

  ! (4) Conserved variables in the two-state (left & right) star regions
  UCstarL(HY_DENS)  = dStarL
  UCstarL(HY_ENER)  = UL(HY_ENER)*numerL/denomL + &
               ((pStar*qStar - totalPresL*velNL))/denomL

  UCstarR(HY_DENS)  = dStarR
  UCstarR(HY_ENER)  = UR(HY_ENER)*numerR/denomR + &
               ((pStar*qStar - totalPresR*velNR))/denomR

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD) /* for MHD */
  UCstarL(HY_MAGX:HY_MAGZ)= Uhll(HY_MAGX:HY_MAGZ)
  UCstarL(HY_ENER) = UCstarL(HY_ENER) &
       -(Bn_hll*dot_product(Uhll(HY_MAGX:HY_MAGZ),Uhll(HY_VELX:HY_VELZ))/Uhll(HY_DENS) &
         -magNL*dot_product(  Vm(HY_MAGX:HY_MAGZ),  Vm(HY_VELX:HY_VELZ)))/denomL

  UCstarR(HY_MAGX:HY_MAGZ)= Uhll(HY_MAGX:HY_MAGZ)
  UCstarR(HY_ENER) = UCstarR(HY_ENER) &
       -(Bn_hll*dot_product(Uhll(HY_MAGX:HY_MAGZ),Uhll(HY_VELX:HY_VELZ))/Uhll(HY_DENS) &
         -magNR*dot_product(  Vp(HY_MAGX:HY_MAGZ),  Vp(HY_VELX:HY_VELZ)))/denomR
#endif


  select case (dir)
  case (DIR_X)
     UCstarL(HY_XMOM) = dStarL*qStar
     UCstarL(HY_YMOM) = UL(HY_YMOM)*numerL/denomL
     UCstarL(HY_ZMOM) = UL(HY_ZMOM)*numerL/denomL

     UCstarR(HY_XMOM) = dStarR*qStar
     UCstarR(HY_YMOM) = UR(HY_YMOM)*numerR/denomR
     UCstarR(HY_ZMOM) = UR(HY_ZMOM)*numerR/denomR

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD) /* for MHD */
     UCstarL(HY_YMOM) = UCstarL(HY_YMOM) - (BxStar*ByStar-Vm(HY_MAGX)*Vm(HY_MAGY))/denomL
     UCstarL(HY_ZMOM) = UCstarL(HY_ZMOM) - (BxStar*BzStar-Vm(HY_MAGX)*Vm(HY_MAGZ))/denomL

     UCstarR(HY_YMOM) = UCstarR(HY_YMOM) - (BxStar*ByStar-Vp(HY_MAGX)*Vp(HY_MAGY))/denomR
     UCstarR(HY_ZMOM) = UCstarR(HY_ZMOM) - (BxStar*BzStar-Vp(HY_MAGX)*Vp(HY_MAGZ))/denomR
#endif

  case (DIR_Y)
     UCstarL(HY_XMOM) = UL(HY_XMOM)*numerL/denomL
     UCstarL(HY_YMOM) = dStarL*qStar
     UCstarL(HY_ZMOM) = UL(HY_ZMOM)*numerL/denomL

     UCstarR(HY_XMOM) = UR(HY_XMOM)*numerR/denomR
     UCstarR(HY_YMOM) = dStarR*qStar
     UCstarR(HY_ZMOM) = UR(HY_ZMOM)*numerR/denomR

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD) /* for MHD */
     UCstarL(HY_XMOM) = UCstarL(HY_XMOM) - (ByStar*BxStar-Vm(HY_MAGY)*Vm(HY_MAGX))/denomL
     UCstarL(HY_ZMOM) = UCstarL(HY_ZMOM) - (ByStar*BzStar-Vm(HY_MAGY)*Vm(HY_MAGZ))/denomL

     UCstarR(HY_XMOM) = UCstarR(HY_XMOM) - (ByStar*BxStar-Vp(HY_MAGY)*Vp(HY_MAGX))/denomR
     UCstarR(HY_ZMOM) = UCstarR(HY_ZMOM) - (ByStar*BzStar-Vp(HY_MAGY)*Vp(HY_MAGZ))/denomR
#endif

  case (DIR_Z)
     UCstarL(HY_XMOM) = UL(HY_XMOM)*numerL/denomL
     UCstarL(HY_YMOM) = UL(HY_YMOM)*numerL/denomL
     UCstarL(HY_ZMOM) = dStarL*qStar

     UCstarR(HY_XMOM) = UR(HY_XMOM)*numerR/denomR
     UCstarR(HY_YMOM) = UR(HY_YMOM)*numerR/denomR
     UCstarR(HY_ZMOM) = dStarR*qStar

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD) /* for MHD */
     UCstarL(HY_XMOM) = UCstarL(HY_XMOM) - (BzStar*BxStar-Vm(HY_MAGZ)*Vm(HY_MAGX))/denomL
     UCstarL(HY_YMOM) = UCstarL(HY_YMOM) - (BzStar*ByStar-Vm(HY_MAGZ)*Vm(HY_MAGY))/denomL

     UCstarR(HY_XMOM) = UCstarR(HY_XMOM) - (BzStar*BxStar-Vp(HY_MAGZ)*Vp(HY_MAGX))/denomR
     UCstarR(HY_YMOM) = UCstarR(HY_YMOM) - (BzStar*ByStar-Vp(HY_MAGZ)*Vp(HY_MAGY))/denomR
#endif

  end select
  ! End of calculating HLLC intermediate states ---------------------------

#ifdef FLASH_UGLM_MHD
  UCstarL(HY_GLMP) = Vm(HY_GLMP)
  UCstarR(HY_GLMP) = Vp(HY_GLMP)
#endif

  ! (5) Finally, calculate HLLC fluxes
  if (SL >= 0.) then
     Fstar(F01DENS_FLUX:hyEndFlux) = FL(F01DENS_FLUX:hyEndFlux)

  elseif ((SL < 0.).and. (qStar >= 0.)) then
     Fstar(F01DENS_FLUX:hyEndFlux-1) = FL(F01DENS_FLUX:hyEndFlux-1) &
          + SL*(UCstarL(HY_DENS:hyEndVar) - UL(HY_DENS:hyEndVar))
     Fstar(HY_P_FLUX) = FL(HY_P_FLUX)

  elseif ((qStar <0.) .and. (SR >= 0.)) then
     Fstar(F01DENS_FLUX:hyEndFlux-1) = FR(F01DENS_FLUX:hyEndFlux-1) &
          + SR*(UCstarR(HY_DENS:hyEndVar) - UR(HY_DENS:hyEndVar))
     Fstar(HY_P_FLUX) = FR(HY_P_FLUX)

  else
     Fstar(F01DENS_FLUX:hyEndFlux) = FR(F01DENS_FLUX:hyEndFlux)

  endif

!!$#ifdef FLASH_UGLM_MHD
!!$  Fstar(F09GLMP_FLUX) = Vm(HY_GLMP)
!!$!print*,Fstar(F09GLMP_FLUX)
!!$#endif


#ifdef FLASH_USM_MHD
  ! Enforce zero for corresponding magnetic field components
  Fstar(F06MAGX_FLUX+dir-1) = 0.
#endif

End Subroutine hy_uhd_HLLC
