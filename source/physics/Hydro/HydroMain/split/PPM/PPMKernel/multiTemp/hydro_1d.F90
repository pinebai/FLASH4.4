!!****if* source/physics/Hydro/HydroMain/split/PPM/PPMKernel/multiTemp/hydro_1d
!!
!! NAME
!!
!!  hydro_1d
!!
!!
!! SYNOPSIS
!!  
!!  call hydro_1d(integer(in) :: blockID,
!!                integer(in) :: numIntCells,
!!                integer(in) :: numCells, 
!!                integer(in) :: guard,
!!                integer(in) :: bcs,
!!                integer(in) :: xyzswp, 
!!                integer(in) :: hy_meshMe, 
!!                real(in)    :: dt, 
!!                real(in)    :: dt_old,                 
!!                integer(in) :: jCell, 
!!                integer(in) :: kCell,                             
!!                integer(in) :: igeom, 
!!                logical(in) :: useGravity,                             
!!                real(in)    :: xbot, 
!!                real(in)    :: xtop,                               
!!                real(in)    :: ybot, 
!!                real(in)    :: ytop, 
!!                real(in)    :: ylft, 
!!                real(in)    :: yrgt,                   
!!                real(in)    :: zlft, 
!!                real(in)    :: zrgt, 
!!                real(in)    :: ugrid(numCells),
!!                real(in)    :: primaryCoord(numCells),
!!                real(in)    :: primaryLeftCoord(numCells),
!!                real(in)    :: primaryRghtCoord(numCells),
!!                real(in)    :: primaryDx(numCells),
!!                real(in)    :: secondCoord(numCells),
!!                real(in)    :: thirdCoord(numCells),
!!                real(in)    :: radialCoord(numCells),
!!                real(in)    :: u(numCells), 
!!                real(in)    :: ut(numCells), 
!!                real(in)    :: utt(numCells), 
!!                real(in)    :: rho(numCells), 
!!                real(in)    :: p(numCells,0:HYDRO_NUM_E_COMPONENTS),
!!                real(in)    :: e(numCells,0:HYDRO_NUM_E_COMPONENTS),
!!                real(in)    :: eintIn(numCells,0:HYDRO_NUM_EINT_COMPONENTS),
!!                real(in)    :: tmp(numCells), 
!!                real(in) :: game(numCells,0:HYDRO_NUM_GAME_COMPONENTS), 
!!                real(inout) :: gamc(numCells),   
!!                real(inout) :: xn(numCells, hy_numXn), 
!!                real(inout) :: utbt(numCells), 
!!                real(inout) :: uttp(numCells),
!!                real(inout) :: utlt(numCells),
!!                real(inout) :: utrt(numCells),               
!!                real(in)    :: shock_multid(numCells),
!!                real(out)   :: dtdx(numCells), 
!!                real(inout) :: areaLeft(numCells),
!!                real(out)   :: area(numCells), 
!!                real(in)    :: cvol(numCells), 
!!                real(inout) :: grav(numCells), 
!!                real(out)   :: ngrav(numCells), 
!!                real(out)   :: fict(numCells),            
!!                real(out)   :: rhoflx(numCells),
!!                real(out)   :: uflx(numCells), 
!!                real(out)   :: pav(numCells, 0:HYDRO_NUM_E_COMPONENTS), 
!!                real(out)   :: utflx(numCells), 
!!                real(out)   :: uttflx(numCells),
!!                real(out)   :: eflx(numCells, 0:HYDRO_NUM_E_COMPONENTS),
!!                real(out)   :: eintflx(numCells, 0:HYDRO_NUM_EINT_COMPONENTS),
!!                real(out)   :: eiaflx(numCells, 0:HYDRO_NUM_EINT_COMPONENTS),
!!                real(out)   :: oneflx(numCells),
!!                real(out)   :: xnflx(numCells,hy_numXn))
!!
!! DESCRIPTION
!!
!!  Compute the 1-d directionally split fluxes through the boundary
!!  of the computational zone using the PPM method.   
!!
!! ARGUMENTS
!!
!!   blockID - my blockID
!!   numIntCells -
!!   numCells - 
!!   guard -  number of guard cells
!!   bcs -
!!   xyzswp - the direction of the sweep
!!   hy_meshMe - the local processor number. 
!!   dt -   current delta t
!!   dt_old - delta t from previous step
!!   jCell -   index that indicates where we are in terms of the
!!             first transversal coordinate, i.e., index for secondCoord
!!   kCell -   index that indicates where we are in terms of the
!!             second transversal coordinate, i.e., index for thirdCoord
!!   igeom -
!!   useGravity - indication if gravitational acceleration should be account for
!!   xbot - 
!!   xtop -                               
!!   ybot - 
!!   ytop - 
!!   ylft - 
!!   yrgt -                   
!!   zlft - 
!!   zrgt -  Values of secondary and third coordinates (at center) for the 1d arrays of 
!!                   cells that are directly neighboring the current 1d array being swept.
!!                   Only used in avisco.
!! 
!!   ugrid - 
!!   primaryCoord -  positions of cell centers of the 1d slice
!!
!!   primaryLeftCoord -  positions of left interfaces
!!
!!   primaryRghtCoord -  positions of right interfaces
!!
!!   primaryDx -  width of cells in the sweep direction
!!
!!   secondCoord -  values of the first transverse coordinate to the sweep direction
!!                   for an x sweep: y coordinates; for a y sweep: x coord; for a z sweep: xcoord
!!
!!   thirdCoord -   the second transverse coordinate to the sweep direction
!!                   for an x sweep: z coordinates; for a y sweep: z coord; for a z sweep: ycoord
!!
!!   radialCoord -  The radial coordinate, no matter whether it's primary, second, or third.
!!                   Actually, this is currently always the IAXIS coordinate.
!!   u - 
!!   ut -
!!   utt -
!!   rho -
!!   p -
!!   e -
!!   tmp - 
!!   game - 
!!   gamc -  
!!   xn - 
!!   utbt - 
!!   uttp -
!!   utlt - 
!!   utrt - 
!!   shock_multid -
!!   dtdx - 
!!   areaLeft -
!!   area - 
!!   cvol -
!!   grav -
!!   ngrav - 
!!   fict -
!!   rhoflx - 
!!   uflx - 
!!   pav -
!!   utflx - 
!!   uttflx - 
!!   eflx - 
!!   eintflx -
!!   eiaflx -       Fluxes for Internal Energies Advected
!!   oneflx -       Fluxes of "Constant One" (i.e., just plain velocity x area)
!!   xnflx -
!!   
!!***

#include "Flash.h"

#ifdef SELE_MSCALAR
#define SELE_IND (SELE_MSCALAR-SPECIES_BEGIN+1)
#else
#define SELE_IND 0
#endif

subroutine hydro_1d (blockID,numIntCells,numCells, guard,bcs,        &
                     xyzswp, hy_meshMe, dt, dt_old,                 &
                     jCell, kCell,                             &
                     igeom, useGravity,                             &
                     xbot, xtop,                               &
                     ybot, ytop, ylft, yrgt,                   &
                     zlft, zrgt, ugrid,                        &
                     primaryCoord ,     &
                     primaryLeftCoord , &
                     primaryRghtCoord , &
                     primaryDx        , &
                     secondCoord      , &
                     thirdCoord       , &
                     radialCoord     , &
                     u, ut, utt, rho, p, e, eintIn, tmp, game, gamc,   &
                     xn, utbt, uttp, utlt, utrt,               &
                     shock_multid,                             &
                     dtdx, areaLeft, area, cvol, grav, ngrav, fict, &
                     rhoflx, uflx, pav, utflx, uttflx,         &
                     eflx, eintflx, eiaflx, oneflx, voFlx, xnflx)


  use Hydro_data, ONLY : hy_dirGeom, hy_movingGrid, hy_useCellAreasForFluxes, &
       hy_epsiln, hy_omg1, hy_omg2
  use Hydro_data, ONLY : hy_cvisc, hy_useCmaAdvection, hy_smallx, hy_smallp
  use Hydro_data, ONLY : hy_numXn, hy_hybridRiemann, hy_updateHydroFluxes
  use Hydro_data, ONLY : hy_ppmEintFluxConstructionMeth,hy_ppmEnerFluxConstructionMeth, &
       hy_ppmEintCFluxConstructionMeth,hy_ppmEnerCFluxConstructionMeth
  use Hydro_data, ONLY : hy_eMassInUAmu
  use Hydro_data, ONLY : hy_3Ttry_Arelated, hy_dbgReconstConsvSele
  use Driver_interface, ONLY : Driver_abortFlash
  use Gravity_interface, ONLY : Gravity_accelOneRow
  use Eos_interface, ONLY : Eos_getAbarZbar
  use hy_ppm_interface, ONLY: hy_ppm_force, hy_ppm_geom, hy_ppm_completeGeomFactors
  use hy_ppm_kernelInterface, ONLY: intrfc, states, rieman, hy_dbgWriteStates
  implicit none
  
#include "constants.h"
#include "Hydro_components.h"

!--arguments-------------------------
  integer, intent(IN) ::  blockID,jCell, kCell, numIntCells,numCells,&
                          xyzswp, hy_meshMe,igeom, guard
  real, intent(IN) :: dt, dt_old

  logical, intent(IN) :: useGravity
  integer,intent(IN),dimension(2,MDIM):: bcs

  real, DIMENSION(numCells), intent(IN) :: rho, u, ut, utt, tmp
  real, DIMENSION(numCells), intent(IN) ::  primaryCoord ,     &
                                            primaryLeftCoord , &
                                            primaryRghtCoord , &
                                            primaryDx        , &
                                            secondCoord      , &
                                            thirdCoord       , &
                                            radialCoord     , &
                                            shock_multid
  real, DIMENSION(numCells,0:HYDRO_NUM_E_COMPONENTS),intent(IN) :: p, e
  real, DIMENSION(numCells,0:HYDRO_NUM_EINT_COMPONENTS),intent(IN) :: eintIn
  real, DIMENSION(numCells), intent(IN)    :: cvol
  real, DIMENSION(numCells), intent(INOUT) :: grav, areaLeft
  real, DIMENSION(numCells), intent(OUT)   :: area
  real, DIMENSION(numCells, hy_numXn), intent(OUT) :: xnflx
  real, DIMENSION(numCells), intent(OUT) :: dtdx, ngrav, fict
  real, DIMENSION(numCells), intent(OUT) :: rhoflx, uflx, &
                                            utflx, uttflx
  real, DIMENSION(numCells), intent(OUT) :: oneflx, voFlx
  real, DIMENSION(numCells,0:HYDRO_NUM_E_COMPONENTS), intent(OUT) :: eflx, pav
  real, DIMENSION(numCells,0:HYDRO_NUM_EINT_COMPONENTS), intent(OUT) :: eintflx
  real, DIMENSION(numCells,0:HYDRO_NUM_EINT_COMPONENTS), intent(OUT) :: eiaflx

  

  real, DIMENSION(numCells, hy_numXn), intent(INOUT) :: xn
  real, intent(IN) :: xbot, xtop, ybot 
  real, intent(IN) :: ytop, ylft, yrgt, zlft, zrgt
  real, intent(IN), DIMENSION(numCells) :: ugrid
  real, DIMENSION(numCells,0:HYDRO_NUM_GAME_COMPONENTS), intent(INout) :: game
  real, intent(INOUT), DIMENSION(numCells) :: gamc, utbt, &
                                              uttp, utlt, utrt


!--locals-------------------------

  real, DIMENSION(numCells, hy_numXn) :: xnav, xnl, xnr

  real, DIMENSION(numCells)    :: rhoav, uav, utav, uttav, &
              &                   rhol, rhor, ul, ur, utl, utr, uttl, uttr, &
              &                   pl, pr, vl, vr, gamcl, gamcr, c, ce, &
              &                   urell, ugrdl, &
              &                   v, dvol, &
              &                   ograv, hgrav
  real, DIMENSION(numCells,0:HYDRO_NUM_EINT_COMPONENTS) :: eint, &
              &                   eintl, eintr, eintAv, gameav
  real, DIMENSION(numCells,0:HYDRO_NUM_E_COMPONENTS) :: etot, &
              &                   etotl, etotr, etotAv
  real, DIMENSION(numCells,1:HYDRO_NUM_E_COMPONENTS) :: pFrac, &
              &                   pFracl, pFracr, pFracAv
  real, DIMENSION(numCells,0:HYDRO_NUM_GAME_COMPONENTS) ::  gamel, gamer


  
  real, DIMENSION(numCells) ::  scrch1, scrchEkin
  real, DIMENSION(numCells) ::  avis
  real, DIMENSION(numCells) ::  xzn, yzn, zzn
  real, DIMENSION(numCells) :: xlold, xrold, dxold, dvold, &
                               alold, aold, arear
  
  integer :: i, n, kk
  integer :: numIntCells4, numIntCells5, numIntCells8
  
  real ::  dtfac, dg
  real :: Ye, ekinElecFrac, smallp
  integer :: eintFluxConstructionMethod, fluxConstructionMethod

!  real :: pres_jump
  integer,dimension(2) :: pos
  
#ifndef RHO_FLUX
  call Driver_abortFlash("[HYDRO_1D] ERROR: rhoflx not defined")
#endif
#ifndef U_FLUX
  call Driver_abortFlash("[HYDRO_1D] ERROR: uflx not defined")
#endif
#ifndef P_FLUX
  call Driver_abortFlash("[HYDRO_1D] ERROR: pflx not defined")
#endif
#ifndef UT_FLUX
  call Driver_abortFlash("[HYDRO_1D] ERROR: utflx not defined")
#endif
#ifndef UTT_FLUX
  call Driver_abortFlash("[HYDRO_1D] ERROR: uttflx not defined")
#endif
#ifndef E_FLUX
  call Driver_abortFlash("[HYDRO_1D] ERROR: eflx not defined")
#endif
#ifndef EINT_FLUX
  call Driver_abortFlash("[HYDRO_1D] ERROR: eintflx not defined")
#endif


                        ! initialize grav arrays to zero
  hgrav(:) = 0.e0       ! need to reconsider whether this
  ngrav(:) = 0.e0       ! is the best place to do this


  numIntCells4 = numIntCells + 4
  numIntCells5 = numIntCells + 5
  numIntCells8 = numIntCells + 8

  !---------------------------
  if (.NOT. hy_useCellAreasForFluxes) then
     call hy_ppm_geom (numIntCells, numCells, jCell, kCell, xyzswp, igeom,  &
          &     areaLeft, arear, area, primaryDx, dvol, &
          &     primaryLeftCoord, primaryRghtCoord,  &
          &     radialCoord, thirdCoord)
  else
     call hy_ppm_completeGeomFactors(numIntCells4, numCells, igeom, primaryDx, &
          radialCoord(jCell), dvol, cvol, &
          area, areaLeft)
  end if

  !---------------------------
  if (useGravity) then
     pos(1)=jCell; pos(2)=kCell
#ifdef GPOL_VAR
#ifndef GPOT_VAR
     call Driver_abortFlash("Shouldn't have gpol defined without gpot")
#endif
#endif

#if defined(GPOT_VAR) && defined(GPOL_VAR) && defined(FLASH_GRAVITY_TIMEDEP)
     ! Gravity implementation defines FLASH_GRAVITY_TIMEDEP -> time-dependent gravity field,
     ! interpolate the acceleration linearly in time (pointwise) - KW
     ograv(:) = 0.e0       ! initialize array to zero
     call Gravity_accelOneRow (pos, xyzswp, blockID, numIntCells8,ograv,GPOL_VAR)
     call Gravity_accelOneRow (pos, xyzswp, blockID, numIntCells8,grav,GPOT_VAR)
     dtfac = dt/dt_old

     do i = 1,numIntCells8
        dg       = dtfac*(grav(i) - ograv(i))
        hgrav(i) = grav(i) + 0.5e0*dg
        ngrav(i) = grav(i) +       dg
     enddo

#else
     ! FLASH_GRAVITY_TIMEDEP not defined -> assume time-independent gravity field.
     ! Also if GPOT_VAR or GPOL_VAR defined -> use current accel without time
     ! interpolation, i.e., handle like time-independent gravity field - KW
     call Gravity_accelOneRow (pos, xyzswp, blockID, numIntCells8,grav)
     hgrav = grav
     ngrav = grav
#endif

  else

     hgrav(:) = 0.e0
     ngrav(:) = 0.e0

  end if

  if (hy_updateHydroFluxes) then

    !---------------------------
    ! For non-cartesian geometries, this call sets fict to have non-zero
    ! values.
    !---------------------------
     call hy_ppm_force(numCells, numIntCells, jCell, kCell, igeom, &
                primaryCoord, radialCoord, thirdCoord, u, ut, utt, fict)


     do i = 2, numIntCells8
        ugrdl(i) = 0.5e0 * (ugrid(i) + ugrid(i-1))
     end do
     ugrdl(1) = ugrdl(2)

     do i = 1, numIntCells8
        v(i)  = 1.e0 / rho(i)
        c(i)  = sqrt (gamc(i) * p(i,0) * rho(i))
        ce(i) = c(i)*v(i)
     end do


    ! Initialization -- grid/geometry factors

     do i = 1, numIntCells8
        dtdx(i) = dt / primaryDx(i)
     enddo

     ! Note - if igeom = 3, 4, or 5, the x coordinate will be a radial
     ! coordinate and the index jCell will always refer to the x direction.

     if (igeom >= PHI_CYL) then
        do i = 1, numIntCells8
           dtdx(i) = dtdx(i) / radialCoord(jCell)
        enddo
     endif
     if (igeom == PHI_SPH) then
        do i = 1, numIntCells8
           dtdx(i) = dtdx(i) / sin(thirdCoord(kCell))
        enddo
     endif

     ! compute the internal energy - first time
     do i = 1,numIntCells8 
        scrchEkin(i) = 0.5e0*(u(i)**2 + ut(i)**2 + utt(i)**2)
        do n = 0,HYDRO_NUM_EINT_COMPONENTS
           select case (n)
           case(0)
              etot(i,n) = scrchEkin(i) !abuse etot array for kinetic energy...
              eint(i,n) = e(i,n) - scrchEkin(i)
           case(1)
#if defined(FLASH_MULTISPECIES) || (defined(SUMY_MSCALAR) && defined(YE_MSCALAR))
              Ye = 0.5  !DEV: The modes using this (ppmEnerCompFlux..>=11) to go away anyway. - KW
#else
              call Eos_getAbarZbar(Ye=Ye,massFrac=xn(i,1:NSPECIES))           !DEV: for now
#endif
              ekinElecFrac = Ye * hy_eMassInUAmu
              etot(i,n) = (1-ekinElecFrac)*scrchEkin(i) !abuse etot array for kinetic energy
!!$              eint(i,n) = e(i,n) - (1-ekinElecFrac)*scrchEkin(i)
              eint(i,n) = eintIn(i,n)
           case(2)
              etot(i,n) = ekinElecFrac * scrchEkin(i) !abuse etot array for kinetic energy
!!$              eint(i,n) = e(i,n) - ekinElecFrac * scrchEkin(i)
              eint(i,n) = eintIn(i,n)
           case default
              etot(i,n) = 0.0 !abuse etot array for kinetic energy
              eint(i,n) = eintIn(i,n)
           end select
           eintFluxConstructionMethod = hy_ppmEintFluxConstructionMeth
           if (n>0) eintFluxConstructionMethod = hy_ppmEintCFluxConstructionMeth
           fluxConstructionMethod = hy_ppmEnerFluxConstructionMeth
           if (n>0) fluxConstructionMethod = hy_ppmEnerCFluxConstructionMeth
           if (fluxConstructionMethod .GE. 0 .AND. fluxConstructionMethod < 10) then
              etot(i,n) = eint(i,n) !abuse etot array for internal energy...
           else if (fluxConstructionMethod .GE. 10 .AND. fluxConstructionMethod < 20) then
              etot(i,n) = e(i,n) !Actually use etot array as its name suggests!
           end if
           select case (eintFluxConstructionMethod)
           case(0,1,4,5)
              if (n /= 3) then  !! ERAD can be zero
                 eint(i,n) = max(eint(i,n),hy_smallp/rho(i)) !primitive (mass-specific)
              end if
           case(2,6)
              eint(i,n) = max(eint(i,n),hy_smallp/rho(i)) * rho(i) !conserved (per-vol)
           case(3,7)
              eint(i,n) = max(eint(i,n),hy_smallp/rho(i)) / max(e(i,0)-scrchEkin(i),hy_smallp/rho(i)) !comp as fraction of eint
           end select
        end do
     end do
     do i = 1,numIntCells8 
        do n = HYDRO_NUM_E_COMPONENTS,0,-1
           !          etot(i,n) = e(i,n)   ! now done above.
           fluxConstructionMethod = hy_ppmEnerFluxConstructionMeth
           if (n>0) fluxConstructionMethod = hy_ppmEnerCFluxConstructionMeth
           if (fluxConstructionMethod .GE. 20 .AND. fluxConstructionMethod < 30) then
              smallp = 0.0
           else
              smallp = hy_smallp
           end if
           select case (fluxConstructionMethod)
           case(0,1,4,5, 11,14,15, 20,21,24,25)
              etot(i,n) = max(etot(i,n),smallp/rho(i)) !primitive (mass-specific)
           case(2,6,     12,16,    22,26)
              etot(i,n) = max(etot(i,n),smallp/rho(i)) * rho(i) !conserved (per-vol)
           case(3,7,     13,17)
              etot(i,n) = max(etot(i,n),smallp/rho(i)) / max(etot(i,0),hy_smallp/rho(i)) !comp as fraction of etot
           case(                   23,27)
              if (etot(i,0) > 0) then
                 etot(i,n) = etot(i,n) / etot(i,0) !comp as fraction of etot
              else
                 etot(i,n) = 0.0
              end if
           end select
        end do
     end do
     do i = 1,numIntCells8 
        do n = 1,HYDRO_NUM_E_COMPONENTS
           pFrac(i,n) = p(i,n)/p(i,0)
        end do
     end do
     do i = 1,numIntCells8 
        do n = 1,HYDRO_NUM_GAME_COMPONENTS
           game(i,n) = p(i,n)/( eintIn(i,n) * rho(i)) + 1
        end do
     end do


! Obtain PPM interpolation coefficients.
     call intrfc(xyzswp,numIntCells,numCells, guard,&
          &      rho,u,ut,utt,p, &
          &      rhol,rhor,ul,ur, &
          &      utl,utr,uttl,uttr, &
          &      pl,pr,vl,vr,gamcl, &
          &      gamcr,game, &
          &      gamel,gamer,gamc,hgrav,&
                 eint,eintl,eintr,etot,etotl,etotr,pFrac,pFracl,pFracr,xn, &
          &      xnl,xnr,v,primaryDx, primaryCoord,tmp)

     ! Determine input states for the Riemann solver.
     call states (numIntCells,numCells,&
                  jCell,igeom, &
                  rho,u,rhol,rhor,ul,ur,&
                  utl,utr,uttl,uttr,p, pl,pr, &
                  gamcl,gamcr,&
                  ugrid,ce,game,gamer,gamc,gamel,&
                  eintl,eintr, &
                  etotl,etotr, &
                  pFracl,pFracr, &
                  xnl,xnr, &
                  dtdx, dt, &
                  primaryCoord, primaryLeftCoord, radialCoord, hgrav, fict)
    ! Solve Riemann problems at zone edges.
     call rieman(numIntCells,numCells, &
                 rhoav,uav,utav,uttav,pav,&
                 urell,ugrdl,dvold,game,gameav,pFracAv,eintAv,etotAv,xnav,primaryCoord)

     gameav(:,HYDRO_NUM_GAME_COMPONENTS+1:min(max(HYDROCOMP_ION,HYDROCOMP_ELE),HYDRO_NUM_EINT_COMPONENTS)) = 5./3.
     gameav(:,max(HYDROCOMP_ION,HYDROCOMP_ELE,HYDRO_NUM_GAME_COMPONENTS)+1:min(HYDROCOMP_RAD,HYDRO_NUM_EINT_COMPONENTS)) = 4./3.

     select case (xyzswp)

     case (SWEEP_X)
        if (bcs(LOW,IAXIS) .eq. REFLECTING) then
           uav  (guard+1) = 0.e0
           urell(guard+1) = 0.e0
        endif
        if (bcs(HIGH,IAXIS) .eq. REFLECTING) then
           uav  (guard+numIntCells+1) = 0.e0
           urell(guard+numIntCells+1) = 0.e0
        endif

     case (SWEEP_Y)
        if (bcs(LOW,JAXIS) .eq. REFLECTING) then
           uav  (guard+1) = 0.e0
           urell(guard+1) = 0.e0
        endif
        if (bcs(HIGH,JAXIS) .eq. REFLECTING) then
           uav  (guard+numIntCells+1) = 0.e0
           urell(guard+numIntCells+1) = 0.e0
        endif

     case (SWEEP_Z)
        if (bcs(LOW,KAXIS) .eq. REFLECTING) then
           uav  (guard+1) = 0.e0
           urell(guard+1) = 0.e0
        endif
        if (bcs(HIGH,KAXIS) .eq. REFLECTING) then
           uav  (guard+numIntCells+1) = 0.e0
           urell(guard+numIntCells+1) = 0.e0
        endif
     end select

     ! -------------------------------------------------------------------
     ! Consistent Multi-fluid Advection (Plewa & Mueller 1999, CMA   Eq. 13)

     ! mass fractions renormalization and optional limiting (not if CMA)

     if ( NSPECIES > 1  ) then

        scrch1(5:numIntCells5) = 0.e0

        do n = 1, NSPECIES

          ! renormalize and limit mass fractions: note that limiting introduces
          ! non-conservation of species

           if ( .not.hy_useCmaAdvection ) then
              do i = 5, numIntCells5
                 xnav(i,n) = max(hy_smallx, min(1.e0, xnav(i,n)))
              end do
           end if

           do i = 5, numIntCells5
              scrch1(i) = scrch1(i) + xnav(i,n)
           end do
        end do

        do i = 5, numIntCells5
           if ( scrch1(i) /= 0.e0 ) scrch1(i) = 1.e0 / scrch1(i)
        end do

        do n = 1, NSPECIES
           do i = 5, numIntCells5
              xnav(i,n) = xnav(i,n) * scrch1(i)
           end do
        end do

     end if

     !--------------------------------------------------------------------------
     ! Save old grid information and move the grid using the previously computed 
     ! grid velocities (for artificial dissipation purposes).

     do i = 1, numIntCells8
!       xlold(i) = primaryLeftCoord(i) !  Use of xlold is commented out below.
!       xrold(i) = primaryRghtCoord(i) !  xrold is never used.
!       dxold(i) = primaryDx(i)  ! dxold is never used.
!       dvold(i) = dvol(i)       ! dvold used to be never used.
        alold(i) = areaLeft(i)
!       aold(i)  = area(i)       ! aold is never used.
     enddo

     if (hy_movingGrid) then !! this is not defined yet
!!$       do i = 2, numIntCells8
!!$          primaryLeftCoord(i)   = xlold(i) + dt * ugrdl(i)
!!$          primaryRghtCoord(i-1) = primaryLeftCoord(i)
!!$       enddo
!!$       
!!$       do i = 2, numIntCells8
!!$          primaryCoord(i)  = 0.5e0 * (primaryRghtCoord(i) + primaryLeftCoord(i))    
!!$          primaryDx(i) = primaryRghtCoord(i) - primaryLeftCoord(i)
!!$       enddo

!!  This 2nd call to geom should not be necessary - the geometry better not
!! have changed since the first call above - unless the grid is moving !
!! The following geom call and the geom implementation may not be correct 
!! for a moving grid anyway.
!!$    call hy_ppm_geom (numIntCells, numCells,jCell, kCell, xyzswp, igeom, &
!!$         &     areaLeft, arear, area, primaryDx, dvol, &
!!$               primaryLeftCoord, primaryRghtCoord,  &
!!$               radialCoord, thirdCoord)
     endif

     !------------------------------------------------------------------------------
     ! Compute unmodified fluxes for each of the conserved quantities.
     !
     !For different ppmEnerFluxConstructionMeth methods, energy fluxes are
     !constructed here from the following rieman outputs:
     !
     !        |                  Rieman results (assuming that urell==uav)
     ! Method |    used for eintflx(*,0)    |    used for eflx(*,0)
     !==============================================================================
     !   0    | rhoav,uav,       pav,gameav | rhoav,uav,       pav,gameav,utav,uttav
     !   1,2* | rhoav,uav,eintAv,pav        | rhoav,uav,eintAv,pav       ,utav,uttav
     !   4    | rhoav,uav,       pav,gameav | rhoav,uav,       pav,gameav,utav,uttav
     !   5,6* | rhoav,uav,eintAv,    gameav | rhoav,uav,eintAv,    gameav,utav,uttav
     !Possible optimization:
     !   0    |       uav,       pav,gameav | rhoav,uav,       pav,gameav,utav,uttav
     !   1,2* | rhoav,uav,eintAv,pav        | rhoav,uav,eintAv,pav       ,utav,uttav
     !   4    |       uav,       pav,gameav | rhoav,uav,       pav,gameav,utav,uttav
     !   5,6* | rhoav,uav,eintAv,    gameav | rhoav,uav,eintAv,    gameav,utav,uttav
     !==============================================================================
     ! * For methods 2 and 6, PPM reconstruction, interpolation, and advection for
     !   eintAv is applied to a conserved internal energy variable (i.e., internal
     !   energy expressed in per-volume form).
     !   For methods 3 and 7, PPM reconstruction, interpolation, and advection for
     !   eintAv is applied to an internal energy fraction.  (This makes only sense for
     !   components i.e., hy_ppmEintCFluxConstructionMeth, not hy_ppmEintFluxConstructionMeth).
     !   Otherwise, the primitive (specific) internal energy is used.
     
     
     do i = 5, numIntCells5
        oneflx(i) =             urell(i)
        rhoflx(i) = rhoav (i) * urell(i)
        voFlx(i)  =             urell(i)/dvold(i)
        uflx(i)   = rhoflx(i) * uav  (i)
        utflx(i)  = rhoflx(i) * utav (i)
        uttflx(i) = rhoflx(i) * uttav(i)
        pav(i,1:)  = pav(i,0) * pFracAv(i,:)
     end do



!***********************************************************************************
!***********************************************************************************

     do n = 0, HYDRO_NUM_EINT_COMPONENTS
        eintFluxConstructionMethod = hy_ppmEintFluxConstructionMeth
        if (n>0) eintFluxConstructionMethod = hy_ppmEintCFluxConstructionMeth
        do i = 5, numIntCells5
           select case (eintFluxConstructionMethod)
           case(0,4)
              scrch1(i) = pav(i,n) / ( rhoav(i) * (gameav(i,n)-1.e0) )
           case(1,5)
              scrch1(i) = eintAv(i,n)
           case(2,6)
              scrch1(i) = eintAv(i,n) / rhoav(i)
           case(3,7)
              scrch1(i) = eintAv(i,n) * eintAv(i,0)
           end select

           ! compute the internal energy flux

          if(hy_3Ttry_Arelated) then
           eiaflx(i,n)  = rhoflx(i) * scrch1(i)
          else
           select case (eintFluxConstructionMethod)
           case(0,4,1,5)
              eiaflx(i,n)  = rhoflx(i) * eintAv(i,n)
           case(2,6)
              eiaflx(i,n)  = oneflx(i) * eintAv(i,n)
           case(3,7)
              eiaflx(i,n)  = rhoflx(i) * eintAv(i,n) * eintAv(i,0)
           end select
          endif
           select case (eintFluxConstructionMethod)
           case(0,1,2,3)
              eintflx(i,n) = rhoflx(i) * scrch1(i) + uav(i) * pav(i,n)
           case(4,5,6,7)
              eintflx(i,n) = rhoflx(i) * scrch1(i) * gameav(i,n) 
           end select
        end do
     end do

       ! add the kinetic energy 
     do n = 0, HYDRO_NUM_E_COMPONENTS
        fluxConstructionMethod = hy_ppmEnerFluxConstructionMeth
        if (n>0) fluxConstructionMethod = hy_ppmEnerCFluxConstructionMeth
        do i = 5, numIntCells5

           if (fluxConstructionMethod .GE. 0 .AND. fluxConstructionMethod < 10) then
!!$       scrch1(i) = scrch1(i) &
!!$            &        + 0.5e0 * (uav(i)**2 + utav(i)**2 + uttav(i)**2)
              scrchEkin(i) = 0.5e0 * (uav(i)**2 + utav(i)**2 + uttav(i)**2)
              select case (n)
              case(0)
                 ! NOP
              case(1)
#if defined(FLASH_MULTISPECIES) || (defined(SUMY_MSCALAR) && defined(YE_MSCALAR))
                 Ye = 0.5  !DEV: The modes using eflx(:,1:2) in hy_ppm_updateSoln to go away anyway. - KW
#else
                 call Eos_getAbarZbar(Ye=Ye,massFrac=xn(i,1:NSPECIES))           !DEV: for now
#endif
                 ekinElecFrac = Ye * hy_eMassInUAmu
                 scrchEkin(i) = (1-ekinElecFrac)*scrchEkin(i)
              case(2)
#if defined(FLASH_MULTISPECIES) || (defined(SUMY_MSCALAR) && defined(YE_MSCALAR))
                 Ye = 0.5  !DEV: The modes using eflx(:,1:2) in hy_ppm_updateSoln to go away anyway. - KW
#else
                 call Eos_getAbarZbar(Ye=Ye,massFrac=xn(i,1:NSPECIES))           !DEV: for now
#endif
                 ekinElecFrac = Ye * hy_eMassInUAmu
                 scrchEkin(i) = ekinElecFrac * scrchEkin(i)
              case default
                 scrchEkin(i) = 0.0
              end select
              ! compute the total energy flux
!!$       eflx(i,0)   = rhoflx(i) * scrch1(i) + uav(i) * pav(i,0)
              eflx(i,n)   = eintflx(i,n) + rhoflx(i) * scrchEkin(i)
           else if (fluxConstructionMethod .GE. 20 .AND. fluxConstructionMethod < 30) then
              select case (fluxConstructionMethod)
              case(20,21,24,25)
                 scrchEkin(i) = etotAv(i,n)
              case(22,26)
                 scrchEkin(i) = etotAv(i,n) / rhoav(i)
              case(23,27)
                 scrchEkin(i) = etotAv(i,n) * etotAv(i,0)
              end select
              eflx(i,n)   = eintflx(i,n) + rhoflx(i) * scrchEkin(i)
           else
              select case (fluxConstructionMethod)
              case(11,14,15)
                 scrch1(i) = etotAv(i,n)
              case(12,16)
                 scrch1(i) = etotAv(i,n) / rhoav(i)
              case(13,17)
                 scrch1(i) = etotAv(i,n) * etotAv(i,0)
              end select

              select case (fluxConstructionMethod)
              case(11,12,13)
                 eflx(i,n) = rhoflx(i) * scrch1(i) + uav(i) * pav(i,n)
              case(14,15,16,17)
                 eflx(i,n) = rhoflx(i) * scrch1(i)
                 eintFluxConstructionMethod = hy_ppmEintFluxConstructionMeth
                 if (n>0) eintFluxConstructionMethod = hy_ppmEintCFluxConstructionMeth
                 select case (eintFluxConstructionMethod)
                 case(14)
                    scrch1(i) = pav(i,n) / ( rhoav(i) * (gameav(i,n)-1.e0) )
                 case(15)
                    scrch1(i) = eintAv(i,n)
                 case(16)
                    scrch1(i) = eintAv(i,n) / rhoav(i)
                 case(17)
                    scrch1(i) = eintAv(i,n) * eintAv(i,0)
                 end select
                 eflx(i,n) = eflx(i,n) + rhoflx(i) * scrch1(i) * ( gameav(i,n) - 1)
              end select
           end if
        enddo
     enddo
!!$999 format(A7,I3,':',9(1P,G16.9,1x))
!!$    print 999,'eintflx',blockID,eintflx(5:numIntCells5)


     ! initialize the artificial viscosity coefficient
     avis = 0.0
     ! recompute the internal energy, over the whole structure
     do n = 0, HYDRO_NUM_EINT_COMPONENTS
        do i = 1,numIntCells8 
           scrchEkin(i) = 0.5e0*(u(i)**2 + ut(i)**2 + utt(i)**2)
           select case (n)
           case(0)
              eint(i,n) = e(i,n) - scrchEkin(i)
           case(1)
              eint(i,n) = eintIn(i,n)
           case(2)
              eint(i,n) = eintIn(i,n)
           case default
              eint(i,n) = eintIn(i,n)
           end select
           if (n .ne. HYDROCOMP_RAD) eint(i,n) = max(eint(i,n),hy_smallp/rho(i))
        enddo
     enddo

     !------------------------------------------------------------------------------
     ! Compute quantities needed for artificial viscosity.  Unfortunately in
     ! multidimensions this requires some direction-dependent code.
     
     
     if (hy_cvisc > 0.e0)  then 

        if (xyzswp == SWEEP_X)  then
           xzn(:) = primaryCoord(:)
           yzn(:) = secondCoord(jCell)
           zzn(:) = thirdCoord(kCell)
        endif
        if (xyzswp == SWEEP_Y)  then
           xzn(:) = secondCoord(jCell)
           yzn(:) = primaryCoord(:)
           zzn(:) = thirdCoord(kCell)
        endif
        if (xyzswp == SWEEP_Z)  then
           xzn(:) = secondCoord(jCell)
           yzn(:) = thirdCoord(kCell)
           zzn(:) = primaryCoord(:)
        endif

        call avisco( jCell, kCell , avis,                              &
             hy_dirGeom, xyzswp, 5, numIntCells5, NDIM,    &
             xtop, xbot, ytop, ybot, ylft, yrgt,           &
             zlft, zrgt,                                       &
             primaryCoord, primaryLeftCoord, xzn, yzn, zzn,     &
             u, uttp, utbt, utrt, utlt, hy_cvisc  )

        ! there should be no flux through a reflecting boundary, so force the 
        ! artificial viscosity there to 0

        if (xyzswp == SWEEP_X) then
           if (bcs(LOW,IAXIS) == &
                REFLECTING) avis(guard+1) = 0.e0
           if (bcs(HIGH,IAXIS) == &
                REFLECTING) avis(guard+numIntCells+1) = 0.e0
        endif

        if (xyzswp == SWEEP_Y) then
           if (bcs(LOW,JAXIS) == &
                REFLECTING) avis(guard+1) = 0.e0
           if (bcs(HIGH,JAXIS) == &
                REFLECTING) avis(guard+numIntCells+1) = 0.e0
        endif

        if (xyzswp == SWEEP_Z) then
           if (bcs(LOW,KAXIS) == &
                REFLECTING) avis(guard+1) = 0.e0
           if (bcs(HIGH,KAXIS) == &
                REFLECTING) avis(guard+numIntCells+1) = 0.e0
        endif

        avis(5:numIntCells5) = avis(5:numIntCells5)*alold(5:numIntCells5)
     endif

     ! odd-even decoupling fix
     !
     ! now that we've called the accurate Riemann solver, loop over all the zones,
     ! and for any that is inside a shock, call the HLLE Riemann solver, and 
     ! replace the fluxes computed above with those.  This fixes the odd-even
     ! decoupling problem (see Quirk 1997)
     
     if (hy_hybridRiemann) then

        do i = 5, numIntCells5

          ! check for the presence of shocks by looking at the artificial viscosity
          ! and the pressure jump

!          pres_jump = abs(p(i) - p(i-1))/min(p(i),p(i-1))

!          if (pres_jump >= dp_sh .AND. avis(i) /= 0.e0) then
           if ((shock_multid(i-1)==1 .OR. shock_multid(i)==1) .AND. avis(i) /= 0.e0) then

              ! the interface between zones i and i-1 is a shock.  Use the HLLE solver
              ! and replace the fluxes computed above

              call riemann_hlle( &
                   rho(i-1), rho(i), &
                   u(i-1), u(i), &
                   ut(i-1), ut(i), &
                   utt(i-1), utt(i), &
                   p(i-1,0), p(i,0), &
                   e(i-1,0), e(i,0), &
                   eint(i-1,0), eint(i,0), &
                   game(i-1,0), game(i,0), &
                   gamc(i-1), gamc(i), &
                   hgrav(i-1), hgrav(i), &
                   dt, &
                   rhoflx(i), uflx(i), utflx(i), uttflx(i), eflx(i,0), eintflx(i,0))
             
           endif
          
        enddo

     endif                      !if (hy_hybridRiemann)

     ! Do the elemental abundance fluxes.

     do n = 1, hy_numXn
        do i = 5, numIntCells5
#if (SELE_IND == 0)
           xnflx(i,n) = xnav(i,n) * rhoflx(i)
#else
           if(hy_dbgReconstConsvSele .AND. n==SELE_IND) then
              xnflx(i,n) = xnav(i,n) * oneflx(i)
           else
              xnflx(i,n) = xnav(i,n) * rhoflx(i)
           end if
#endif
        enddo
     enddo

     !---------------------------------------------------------------
     ! Apply the diffusive fluxes to the unmodified conserved fluxes.

     do i = 5, numIntCells5
        oneflx(i) = oneflx(i) * alold(i)
        voFlx(i)  = voFlx(i)  * alold(i)
        rhoflx(i) = rhoflx(i) * alold(i)
        uflx  (i) = uflx  (i) * alold(i)
        utflx (i) = utflx (i) * alold(i)
        uttflx(i) = uttflx(i) * alold(i)
        eflx(i,:) = eflx(i,:) * alold(i)

        eintflx(i,:) = eintflx(i,:)*alold(i)
        eiaflx (i,:) = eiaflx (i,:)*alold(i)

        rhoflx(i)  = rhoflx(i)  + avis(i) * (rho(i-1)           - rho(i))  
        uflx  (i)  = uflx  (i)  + avis(i) * (rho(i-1)*u(i-1)    - rho(i)*u(i))
        utflx (i)  = utflx (i)  + avis(i) * (rho(i-1)*ut(i-1)   - rho(i)*ut(i))
        uttflx(i)  = uttflx(i)  + avis(i) * (rho(i-1)*utt(i-1)  - rho(i)*utt(i))
     enddo

     do n = 0, HYDRO_NUM_E_COMPONENTS
        do i = 5, numIntCells5
           eflx(i,n)  = eflx(i,n)  + avis(i) * (rho(i-1)*e(i-1,n)  - rho(i)*e(i,n))
        end do
     end do

     do n = 0, HYDRO_NUM_EINT_COMPONENTS
        do i = 5, numIntCells5
           eintflx(i,n) = eintflx(i,n) + avis(i) * (rho(i-1)*eint(i-1,n) - rho(i)*eint(i,n))
           eiaflx (i,n) = eiaflx (i,n) + avis(i) * (rho(i-1)*eint(i-1,n) - rho(i)*eint(i,n))
        end do
     end do

     do n = 1, hy_numXn
        do i = 5, numIntCells5
           xnflx(i,n) = xnflx(i,n) * alold(i)
           xnflx(i,n) = xnflx(i,n) + &
                avis(i) * (rho(i-1)*xn(i-1,n) - rho(i)*xn(i,n))
        enddo
     enddo

#ifdef HYDRO_DEBUG_WRITESTATES
    if (hy_meshMe == 0) then
     call hy_dbgWriteStates('rho',rho,rhol,rhor, &
          rhoav,rhoflx,&
          rho,rhol,rhor, rhoav,&
          urell,dt,&
          numIntCells=numIntCells, numCells=numCells, nguard=guard,&
          xCoords=primaryCoord,xL=primaryLeftCoord,xR=primaryRghtCoord,&
          blockID=blockID)
     call hy_dbgWriteStates('etot',etot(:,0),etotl(:,0),etotr(:,0), &
          etotav(:,0),eflx(:,0),&
          rho,rhol,rhor, rhoav,&
          urell,dt,&
          numIntCells=numIntCells, numCells=numCells, nguard=guard,&
          xCoords=primaryCoord,xL=primaryLeftCoord,xR=primaryRghtCoord,&
          blockID=blockID)
     call hy_dbgWriteStates('eint',eint(:,0),eintl(:,0),eintr(:,0), &
          eintav(:,0),eintflx(:,0),&
          rho,rhol,rhor, rhoav,&
          urell,dt,&
          numIntCells=numIntCells, numCells=numCells, nguard=guard,&
          xCoords=primaryCoord,xL=primaryLeftCoord,xR=primaryRghtCoord,&
          blockID=blockID)
     call hy_dbgWriteStates('eion',eint(:,1),eintl(:,1),eintr(:,1), &
          eintav(:,1),eintflx(:,1),&
          rho,rhol,rhor, rhoav,&
          urell,dt,&
          numIntCells=numIntCells, numCells=numCells, nguard=guard,&
          xCoords=primaryCoord,xL=primaryLeftCoord,xR=primaryRghtCoord,&
          blockID=blockID)
     call hy_dbgWriteStates('eele',eint(:,2),eintl(:,2),eintr(:,2), &
          eintav(:,2),eintflx(:,2),&
          rho,rhol,rhor, rhoav,&
          urell,dt,&
          numIntCells=numIntCells, numCells=numCells, nguard=guard,&
          xCoords=primaryCoord,xL=primaryLeftCoord,xR=primaryRghtCoord,&
          blockID=blockID)
#if SELE_IND > 0
     if (hy_dbgReconstConsvSele) then
        call hy_dbgWriteStates('sele',xn(:,SELE_IND)*rho(:),&
          xnl(:,SELE_IND),xnr(:,SELE_IND), &
          xnav(:,SELE_IND),xnflx(:,SELE_IND),&
          rho,rhol,rhor, rhoav,&
          urell,dt,&
          numIntCells=numIntCells, numCells=numCells, nguard=guard,&
          xCoords=primaryCoord,xL=primaryLeftCoord,xR=primaryRghtCoord,&
          blockID=blockID)
     else
        call hy_dbgWriteStates('sele',xn(:,SELE_IND),&
          xnl(:,SELE_IND),xnr(:,SELE_IND), &
          xnav(:,SELE_IND),xnflx(:,SELE_IND),&
          rho,rhol,rhor, rhoav,&
          urell,dt,&
          numIntCells=numIntCells, numCells=numCells, nguard=guard,&
          xCoords=primaryCoord,xL=primaryLeftCoord,xR=primaryRghtCoord,&
          blockID=blockID)
     end if
#endif
    end if
#endif

  else
     fict = 0.0
  end if   !if (hy_updateHydroFluxes)
  !------------------------------------------------------------------------------
  ! Store dt/dx, geometry factors, and the modified fluxes in 'global' arrays
  ! for use in updating the solution after all of the 1D sweeps in this direction
  ! are done. Note that dt/dx is not constant in non-Cartesian geometries, so it
  ! is saved for each zone in the tempDtDx() array.

  do i = 5, numIntCells4
     if (dvol(i) == 0.) then
        print *,' ERROR: dvol == 0 '
     end if
     dtdx(i) = dt/dvol(i)
  enddo

!===================================================================
    
  return
end subroutine hydro_1d
