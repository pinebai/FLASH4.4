!!****if* source/physics/Hydro/HydroMain/split/PPM/multiTemp/hy_ppm_updateSoln
!! NAME
!!
!!  hy_ppm_updateSoln
!!
!!
!! SYNOPSIS
!!
!!  hy_ppm_updateSoln(integer, intent(IN)  :: rangeSwitch, 
!!                integer, intent(IN)  :: xyzswp, 
!!                real, intent(IN)     :: dt,           
!!                integer, intent(IN)  :: blkLimits(HIGH,MDIM),
!!                integer, intent(IN)  :: blkLimitsGC(HIGH,MDIM),
!!                integer, intent(IN)  :: numCells,  
!!                real, intent(IN)     :: tempArea(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
!!                                                 blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
!!                                                 blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),
!!                real, intent(IN)     :: tempGrav1d_o(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
!!                                                     blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
!!                                                     blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),
!!                real, intent(IN)     :: tempGrav1d(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
!!                                                   blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
!!                                                   blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),
!!                real, intent(IN)     :: tempDtDx(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
!!                                                 blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
!!                                                 blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),
!!                real, intent(IN)     :: tempFict(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
!!                                                 blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
!!                                                 blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),
!!                real, intent(IN)     :: tempFlx(NFLUXES, &
!!                                                 blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
!!                                                 blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
!!                                                 blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),
!!                real, pointer,       :: solnData(:,:,:,:)
!!
!! DESCRIPTION
!!
!!  Update the cell average quantities in time based on the fluxes through the 
!!  boundaries.  This is a general update routine for conservation laws in a 
!!  directionally split, finite-volume method.  Given the time-averaged, 
!!  cell-edge fluxes (tempFlx,...), update the solution vector to the new
!!  time in the direction specified by xyzswp.
!!
!!  rangeSwitch specifies which cells to update.  When performing flux 
!!  conservation, it is sometimes useful to delay the update at the boundaries
!!  until after the corrected fluxes have been computed.  rangeSwitch lets
!!  you update either all the cells, all the cells away from the boundary,
!!  or only the cells abutting a boundary.  
!!
!! ARGUMENTS
!!
!!  rangeSwitch --   If UPDATE_INTERIOR, update only "interior" cells. That means here,
!!                    all cells except the first and last layer in the sweep direction.
!!                   If UPDATE_BOUND, update only the first and the last layer of cells
!!                    in the sweep direction.
!!                   Otherwise update all cells.
!!                   Guard cells are never updated.
!!
!!  xyzswp --        The direction of the current sweep, one of SWEEP_X, SWEEP_Y, SWEEP_Z
!!
!!  dt --            The timestep to advance the solution by
!!  blkLimits  -- Index limits of the block interior
!!  blkLimitsGC -- Index limits of the block including guardcell
!!  numCells  -- number of cells in one dimension
!!
!!  tempArea --      Temp data from hydro_1d
!!  tempGrav1d_o -
!!  tempGrav1d -
!!  tempDtDx -       The timestep divided by the cell volume (with 
!!                   geometrical factors)
!!  tempFict -       Geometry related forces, e.g., centrifugal; 
!!                   used to update velocities 
!!
!!  tempFlx --       The fluxes through the boundary
!!
!!  solnData --      Pointer to a block of data
!!
!!
!!***

! solnData depends on the ordering on unk
!!REORDER(4): solnData, tempFlx

#include "Flash.h"

#define CIP_
#ifdef DEBUG_ALL
#define DEBUG_HYDRO
#endif

subroutine hy_ppm_updateSoln(rangeSwitch,                        &
                         xyzswp, dt,                          &
                         blkLimits,blkLimitsGC,numCells,               &
                         tempArea, tempGrav1d_o, tempGrav1d,  &
                         tempDtDx, tempFict,                  &
                         tempFlx,  solnData )
!==============================================================================

  use Hydro_interface, ONLY: Hydro_recalibrateEints
  use Hydro_data, ONLY: hy_numXn
  use Hydro_data, ONLY: hy_smlrho, hy_smallp, hy_eintSwitch, hy_useCmaAdvection,&
       hy_eint1Switch,hy_eint2Switch,hy_eint3Switch
!!$  use Hydro_data, ONLY : hy_eMass, hy_pMass !now unused
  use Hydro_data, ONLY : hy_eMassInUAmu
  use Hydro_data, ONLY : hy_3Ttry_B, hy_3Ttry_D, hy_3Ttry_E, hy_3Ttry_F, hy_3Ttry_G, &
       hy_3Ttry_B_rad, hy_3Ttry_useShockDetect
#ifdef FLASH_MULTISPECIES
  use Eos_interface, ONLY : Eos_getAbarZbar
#else
#if defined(SUMY_MSCALAR) && defined(YE_MSCALAR)
#else
  use Eos_interface, ONLY : Eos_getAbarZbar
#endif
#endif

  implicit none  
#include "constants.h"
#include "PPM.h"

  integer, intent(IN) :: rangeSwitch
  integer, intent(IN) :: xyzswp
  real,    intent(IN) :: dt
  integer,intent(IN) :: numCells  
  integer, intent(IN),dimension(2,MDIM)::blkLimitsGC,blkLimits
#ifdef FIXEDBLOCKSIZE
  real, intent(IN), DIMENSION(GRID_ILO_GC:GRID_IHI_GC, &
                              GRID_JLO_GC:GRID_JHI_GC, &
                              GRID_KLO_GC:GRID_KHI_GC  ) :: &
                                                  tempArea, tempGrav1d_o, &
                                                  tempGrav1d, &
                                                  tempDtDx, tempFict
  real, intent(IN), DIMENSION(NFLUXES,GRID_ILO_GC:GRID_IHI_GC, &
                        GRID_JLO_GC:GRID_JHI_GC, &
                        GRID_KLO_GC:GRID_KHI_GC) :: tempFlx
#else
  real, intent(IN), DIMENSION(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                              blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                              blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)  ) :: &
                                                  tempArea, tempGrav1d_o, &
                                                  tempGrav1d, &
                                                  tempDtDx, tempFict
  real, intent(IN), DIMENSION(NFLUXES,blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                        blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                        blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) :: tempFlx

#endif

  real, pointer ::  solnData(:,:,:,:) 

  integer :: i, j, k, kk, n
  integer :: imin, imax, iskip, jmin, jmax, jskip, kmin, kmax, kskip
  real    :: xnflx2(hy_numXn), xnflx1(hy_numXn), dtdx(numCells)

  real    :: rhoflx1, rhoflx2, uflx1, uflx2, pav1, pav2, utflx1, &
             utflx2, uttflx1, uttflx2, eflx1, eflx2, eintflx1, eintflx2
  real    :: e1flx1,e1flx2, e2flx1,e2flx2, e3flx1,e3flx2
  real    :: eint1flx1,eint1flx2, eint2flx1,eint2flx2, eint3flx1,eint3flx2
  real    :: eiaflx1, eiaflx2
  real    :: oneflx1, oneflx2
  real    :: voFlx(2)
!!$  real    :: transVolV(1:2)
  real,dimension(0:2) :: prevMassV,prevVolV,newVolV,kappa
  real,dimension(0:3,0:2) :: prevEcompV
  real,dimension(0:3,0:2) :: newEcompV
  real,dimension(0:3) :: gamComp,newEinternalV
  integer :: side,piece
  real    :: eia1flx1,eia1flx2, eia2flx1,eia2flx2, eia3flx1,eia3flx2
  logical :: inShock, forceRageLike
  real    :: p1av1,p1av2, p2av1,p2av2, p3av1,p3av2
  real    :: Ye, ekinElecFrac
  real    :: eintIncreaseV, einternalAdvectedV, eintIncreaseAboveAdvectedV
  real    :: eEleAdvectedV, eIonAdvectedV, eRadAdvectedV, eIonEleAdvectedV
  real    :: eEleAdvectedPlusPdVV, eIonAdvectedPlusPdVV, eRadAdvectedPlusPdVV
  real    :: eIonEleAdvectedPlusPdVV, zeroDummy
  real    :: rescaledEEleAdvectedPlusPdVV,rescaledEIonAdvectedPlusPdVV,rescaledERadAdvectedPlusPdVV
  real    :: ekinIonFrac
  real    :: eionIncreaseAboveAdvectedV, eionTargetIncreaseV
  real    :: eionAdjustment, eeleAdjustment

  real    :: aold_t, grav1d_o, grav1d, fict1d
  real    :: aux1, ekin, einternal, etot
  real    :: rho_o, velx_o, vely_o, velz_o, inv_new_dens
  real    :: gamE!, gamERad
  real    :: PiP,PeP,PrP


#ifdef CIP
  integer, save :: itrcr
  real, save    :: pi
  real          :: rho_av, trcr_av
#endif

!===============================================================================
  gamE = 5.0/3.0
  gamComp(0) = gamE
  gamComp(1) = 5.0/3.0
  gamComp(2) = 5.0/3.0
  gamComp(3) = 4.0/3.0

! update in x direction
  if (xyzswp == SWEEP_X) then
     
     select case (rangeSwitch)
     case (UPDATE_INTERIOR)
        imin  = blkLimits(LOW,IAXIS)+ 1
        imax  = blkLimits(HIGH,IAXIS) - 1
        iskip = 1
     case (UPDATE_BOUND)
        imin  = blkLimits(LOW,IAXIS)
        imax  = blkLimits(HIGH,IAXIS)
        iskip = blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)
     case default
        imin  = blkLimits(LOW,IAXIS)
        imax  = blkLimits(HIGH,IAXIS)
        iskip = 1
     end select
     
#ifdef DEBUG_HYDR
     print*,'the sweep direction is',xyzswp      ! within DEBUG
     print*,'the blkLimits is',blkLimits         ! within DEBUG
     print*,'the blkLimitsGC is',blkLimitsGC     ! within DEBUG
     print*,'the update mode is',rangeSwitch,' numCells',numCells  ! within DEBUG
     print*,'imin etc',imin,imax,iskip           ! within DEBUG
#endif

     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do i = imin, imax, iskip
              dtdx(i)  = tempDtDx(i,j,k)
              
              rhoflx1  = tempFlx(RHO_FLUX,i,j,k)
              rhoflx2  = tempFlx(RHO_FLUX,i+1,j,k)

              uflx1    = tempFlx(U_FLUX,i,j,k)
              uflx2    = tempFlx(U_FLUX,i+1,j,k)

              pav1     = tempFlx(P_FLUX,i,j,k)
              pav2     = tempFlx(P_FLUX,i+1,j,k)

              p1av1     = tempFlx(PION_FLUX,i,j,k)
              p1av2     = tempFlx(PION_FLUX,i+1,j,k)
              p2av1     = tempFlx(PELE_FLUX,i,j,k)
              p2av2     = tempFlx(PELE_FLUX,i+1,j,k)
              p3av1     = tempFlx(PRAD_FLUX,i,j,k)
              p3av2     = tempFlx(PRAD_FLUX,i+1,j,k)

              utflx1   = tempFlx(UT_FLUX,i,j,k)
              utflx2   = tempFlx(UT_FLUX,i+1,j,k)

              uttflx1  = tempFlx(UTT_FLUX,i,j,k)
              uttflx2  = tempFlx(UTT_FLUX,i+1,j,k)

              eflx1    = tempFlx(E_FLUX,i,j,k)
              eflx2    = tempFlx(E_FLUX,i+1,j,k)

              e1flx1    = tempFlx(E1_FLUX,i,j,k)
              e1flx2    = tempFlx(E1_FLUX,i+1,j,k)
              e2flx1    = tempFlx(E2_FLUX,i,j,k)
              e2flx2    = tempFlx(E2_FLUX,i+1,j,k)
              e3flx1    = tempFlx(E3_FLUX,i,j,k)
              e3flx2    = tempFlx(E3_FLUX,i+1,j,k)

              eintflx1 = tempFlx(EINT_FLUX,i,j,k)
              eintflx2 = tempFlx(EINT_FLUX,i+1,j,k)

              eint1flx1 = tempFlx(EION_FLUX,i,j,k)
              eint1flx2 = tempFlx(EION_FLUX,i+1,j,k)
              eint2flx1 = tempFlx(EELE_FLUX,i,j,k)
              eint2flx2 = tempFlx(EELE_FLUX,i+1,j,k)
              eint3flx1 = tempFlx(ERAD_FLUX,i,j,k)
              eint3flx2 = tempFlx(ERAD_FLUX,i+1,j,k)

              eiaflx1 = tempFlx(EIA_FLUX,i,j,k)
              eiaflx2 = tempFlx(EIA_FLUX,i+1,j,k)

              eia1flx1 = tempFlx(EI1A_FLUX,i,j,k)
              eia1flx2 = tempFlx(EI1A_FLUX,i+1,j,k)
              eia2flx1 = tempFlx(EI2A_FLUX,i,j,k)
              eia2flx2 = tempFlx(EI2A_FLUX,i+1,j,k)
              eia3flx1 = tempFlx(EI3A_FLUX,i,j,k)
              eia3flx2 = tempFlx(EI3A_FLUX,i+1,j,k)

              oneflx1 = tempFlx(ONE_FLUX,i,j,k)
              oneflx2 = tempFlx(ONE_FLUX,i+1,j,k)

#ifndef VOLD_FLUX
#define VOLD_FLUX ONE_FLUX
#endif
#ifdef VOLD_FLUX
              voFlx(1) =  tempFlx(VOLD_FLUX,i,j,k)
              voFlx(2) = -tempFlx(VOLD_FLUX,i+1,j,k)
#endif

              do kk = 1,hy_numXn
                 xnflx1(kk) = tempFlx(SPECIES_FLUX_BEGIN+kk-1,i,j,k)
                 xnflx2(kk) = tempFlx(SPECIES_FLUX_BEGIN+kk-1,i+1,j,k)
              end do
#ifdef CIP_DEBUG
              ! recover proper flux for CIP

              if ( itrcr > 0 .and. iskip /= 1 ) then
                 write(*,*) 'flux update soln ',itrcr,itrcr-SPECIES_BEGIN+1,hy_numXn
                 do n = itrcr,itrcr
                    kk = n - SPECIES_BEGIN + 1
                    if ( rhoflx1 /= 0.e0 ) then
                       rho_av  = uflx1/rhoflx1
                       trcr_av = xnflx1(kk)/rho_av
                       trcr_av = atan(trcr_av)/(0.9999d0*pi) + 0.5e0
                       xnflx1(kk) = rhoflx1 * trcr_av
                    else
                       xnflx1(kk) = 0.e0
                    end if
                    if ( rhoflx2 /= 0.e0 ) then
                       rho_av  = uflx2/rhoflx2
                       trcr_av = xnflx2(kk)/rho_av
                       trcr_av = atan(trcr_av)/(0.9999d0*pi) + 0.5e0
                       xnflx2(kk) = rhoflx2 * trcr_av
                    else
                       xnflx2(kk) = 0.e0
                    end if
                 end do
              end if
#endif

              aold_t    = tempArea(i,j,k)
              grav1d_o  = tempGrav1d_o(i,j,k)
              grav1d    = tempGrav1d(i,j,k)
              fict1d    = tempFict(i,j,k)
              
!!$              transVolV(1) =  dtdx(i) * oneflx1
!!$              transVolV(2) = -dtdx(i) * oneflx2 

              rho_o     = solnData(DENS_VAR,i,j,k)
              velx_o    = solnData(VELX_VAR,i,j,k)
              
              if ( hy_useCmaAdvection .and. NSPECIES > 1  ) then
                 ! update the partial mass densities and passive scalars * density

                 prevMassV(1) = 0.0; prevMassV(2) = 0.0
                 do n = 1, hy_numXn
                    solnData(SPECIES_BEGIN-1+n,i,j,k) =               &
                              rho_o*solnData(SPECIES_BEGIN-1+n,i,j,k) &
                            -dtdx(i)*(xnflx2(n) - xnflx1(n))
                    prevMassV(1) = prevMassV(1) + dtdx(i) * xnflx1(n)
                    prevMassV(2) = prevMassV(2) - dtdx(i) * xnflx2(n)
                 end do

                 ! update the total mass density
                 
                 solnData(DENS_VAR,i,j,k) = 0.e0

                 do n = 1, NSPECIES
                    solnData(DENS_VAR,i,j,k) =                        &
                    solnData(DENS_VAR,i,j,k) +                        &
                    solnData(SPECIES_BEGIN-1+n,i,j,k)
                 end do
                 ! recover partial densities and passive scalars * density

                 inv_new_dens = 1.e0/solnData(DENS_VAR,i,j,k)

                 do n = 1, hy_numXn
                    solnData(SPECIES_BEGIN-1+n,i,j,k) =               &
                    solnData(SPECIES_BEGIN-1+n,i,j,k) * inv_new_dens
                 end do

              else
                 ! update the density               
                 solnData(DENS_VAR,i,j,k) = solnData(DENS_VAR,i,j,k) - &
                                       dtdx(i) * (rhoflx2-rhoflx1)

                 prevMassV(1) =  dtdx(i) * rhoflx1
                 prevMassV(2) = -dtdx(i) * rhoflx2


                 ! update the mass fractions and passive scalars
                 do n = 1, hy_numXn

                    solnData(SPECIES_BEGIN-1+n,i,j,k) =               &
                           ( rho_o*solnData(SPECIES_BEGIN-1+n,i,j,k)  &
                            -dtdx(i)*(xnflx2(n) - xnflx1(n))         &
                           )/solnData(DENS_VAR,i,j,k)
                 end do

              end if

              ! limit the density 
              solnData(DENS_VAR,i,j,k) =                              &
                   max(hy_smlrho, solnData(DENS_VAR,i,j,k))
              prevMassV(0) = rho_o  !!! solnData(DENS_VAR,i,j,k)
             
              inv_new_dens = 1.e0/solnData(DENS_VAR,i,j,k)

!#ifdef CIP
!              if ( itrcr > 0 ) then
!                 ! CIP update and limit tracer
!                 do n = itrcr,itrcr
!                    solnData(n,i,j,k) = atan(solnData(n,i,j,k))/(0.9999d0*pi) + 0.5e0
!                    solnData(n,i,j,k) = max(0.e0, min(1.e0, solnData(n,i,j,k) ))
!                 end do
!              end if
!#endif

              ! update the velocities               
              aux1 = -dtdx(i)*( (uflx2 - uflx1)                      &
                               +aold_t*(pav2 - pav1))                &
               +0.5e0*dt                                             &
               *( (solnData(DENS_VAR,i,j,k) + rho_o)*fict1d           &
                 +(rho_o*grav1d_o + solnData(DENS_VAR,i,j,k)*grav1d)  &
                )

              solnData(VELX_VAR,i,j,k) =                              &
                         ( rho_o * solnData(VELX_VAR,i,j,k)           &
                          +aux1                                      &
                         )*inv_new_dens
              
              solnData(VELY_VAR,i,j,k) =                              &
                         ( rho_o * solnData(VELY_VAR,i,j,k)           &
                          -dtdx(i) * (utflx2 - utflx1)               &
                         )*inv_new_dens
              
              solnData(VELZ_VAR,i,j,k) =                              &
                         ( rho_o * solnData(VELZ_VAR,i,j,k)           &
                          -dtdx(i) * (uttflx2 - uttflx1)             &
                         )*inv_new_dens

              ! update the total energy
              aux1 = - dtdx(i) * (eflx2 - eflx1)                             &
                   + dt*0.5e00                                               &
                   *( rho_o*velx_o*grav1d_o                                  &
                     +solnData(DENS_VAR,i,j,k)*solnData(VELX_VAR,i,j,k)*grav1d &
                    )

              etot = (rho_o * solnData(ENER_VAR,i,j,k) + aux1)*inv_new_dens
              
#ifdef EINT_VAR
              ! get the internal energy
              einternal = solnData(EINT_VAR,i,j,k)
              prevEcompV(0,0) = rho_o * solnData(EINT_VAR,i,j,k)

              ! update internal energy
              aux1 = -dtdx(i) * ( (eintflx2 - eintflx1)               &
                     -0.5e0*(velx_o + solnData(VELX_VAR,i,j,k))        &
                     *aold_t * (pav2 - pav1) )
              
              einternal = (rho_o * einternal + aux1)*inv_new_dens
              einternal = max(einternal, hy_smallp*inv_new_dens)

! compute the new kinetic energy       
              ekin = 0.5e0 * ( solnData(VELX_VAR,i,j,k)**2             &
                              +solnData(VELY_VAR,i,j,k)**2             &
                              +solnData(VELZ_VAR,i,j,k)**2)

! test whether we should use the internal energy from the evolution
              if (einternal .LT. hy_eintSwitch*ekin) then
#ifdef DEBUG_MAR201
                 print*,'Updating etot at',i,' from',etot,' to',einternal + ekin
#endif
                 etot = einternal + ekin
              else
                 einternal = etot - ekin
              endif

!!$              gamE = 0.5*(pav2+pav1) * inv_new_dens / einternal
!!$              gamE = 0.5* (game + solnData(PRES_VAR,i,j,k) / (rho_o * solnData(EINT_VAR,i,j,k)))
              gamE = 5.0/3.0
              gamComp(0) = gamE
              gamComp(1) = 5.0/3.0
              gamComp(2) = 5.0/3.0
              gamComp(3) = 4.0/3.0

#ifdef DBGS_VAR
              solndata(DBGS_VAR,i,j,k) = 0.5 * &
                   ( solndata(DBGS_VAR,i,j,k)  + tempFlx(SHOK_FLUX,i,j,k)+tempFlx(SHOK_FLUX,i+1,j,k) )
#endif
              inShock = ( tempFlx(SHOK_FLUX,i,j,k)+tempFlx(SHOK_FLUX,i+1,j,k) > 0 )
              if (hy_3Ttry_useShockDetect) then
                 forceRageLike = (.NOT. inShock)
              else
                 forceRageLike = .FALSE.
              end if

              eintIncreaseV = solnData(DENS_VAR,i,j,k)*einternal - rho_o*solnData(EINT_VAR,i,j,k)
              ! update internal energy advected
              aux1 = -dtdx(i) * (eiaflx2 - eiaflx1)
              einternalAdvectedV = aux1
              eintIncreaseAboveAdvectedV =  eintIncreaseV - einternalAdvectedV
              if (inShock) then
#ifdef DEBUG_XHYDRO
999              format(1x,I2,' inShock Comb:',3(1x,1PG22.15))
                 print 999,i,eintIncreaseV,einternalAdvectedV,eintIncreaseAboveAdvectedV
#endif
                 if (eintIncreaseAboveAdvectedV .LE. 0.0) then
                    !print*,i,' not in shock after all.'
                    inShock = .FALSE.
!#ifdef DEBUG_XHYDRO
                 else
                    !print*,i,' still seems we are in shock...'
!#endif
                 end if
              end if

              solnData(ENER_VAR,i,j,k) = etot
              solnData(EINT_VAR,i,j,k) = einternal
#else
              
              solnData(ENER_VAR,i,j,k) = etot
#endif
              
800           format(a,6(1PG23.16))
!!$             if(hy_3Ttry_B==0 .OR. hy_3Ttry_B==1 .OR. (hy_3Ttry_D.GE.2.0 .AND. hy_3Ttry_E==1)) then
              aux1 = -dtdx(i) * (eia3flx2 - eia3flx1)
              eRadAdvectedV = aux1
              aux1 = -dtdx(i) * (eia2flx2 - eia2flx1)
              eEleAdvectedV = aux1
              aux1 = -dtdx(i) * (eia1flx2 - eia1flx1)
              eIonAdvectedV = aux1
!!$             else if(hy_3Ttry_B_rad==0 .OR. hy_3Ttry_B_rad==1) then
!!$              aux1 = -dtdx(i) * (eia3flx2 - eia3flx1)
!!$              eRadAdvectedV = aux1
!!$             end if


             if(hy_3Ttry_B==3 .OR. hy_3Ttry_B_rad==3) then
                prevVolV(0) = 1.0
                prevVolV(1) = dtdx(i) * voFlx(1)
                prevVolV(2) = dtdx(i) * voFlx(2)
                ! Note prevEcompV(0,0) was initialized above from solnData(EINT_VAR,i,j,k), before that changed.
                prevEcompV(1,0) = rho_o * solnData(EION_VAR,i,j,k)
                prevEcompV(2,0) = rho_o * solnData(EELE_VAR,i,j,k)
                prevEcompV(3,0) = rho_o * solnData(ERAD_VAR,i,j,k)
                do side = 1,2      !first left then right interface
                   select case (side)
                   case(1)
                      prevEcompV(0,side) = dtdx(i) * eiaflx1
                      prevEcompV(1,side) = dtdx(i) * eia1flx1
                      prevEcompV(2,side) = dtdx(i) * eia2flx1
                      prevEcompV(3,side) = dtdx(i) * eia3flx1
                   case(2)
                      prevEcompV(0,side) = - dtdx(i) * eiaflx2
                      prevEcompV(1,side) = - dtdx(i) * eia1flx2
                      prevEcompV(2,side) = - dtdx(i) * eia2flx2
                      prevEcompV(3,side) = - dtdx(i) * eia3flx2
                   end select
                   if (voFlx(side) .GE. 0.0) then ! Stuff is coming in on this side
                      ! ... keep it. 
                   else                           ! Stuff is going out on this side
                      prevMassV(0) = prevMassV(0) + prevMassV(side)
                      prevVolV(0)  = prevVolV(0)  + prevVolV(side)
                      prevEcompV(0:3,0)  = prevEcompV(0:3,0)  + prevEcompV(0:3,side)
                      prevMassV(side) = 0.0
                      prevVolV(side) = 0.0
!!$                    prevEcompV(0:3,side) = 0.0
                   end if
                end do

                newEinternalV(0:3) = 0.0
                do piece = 0,2      !all fluid pieces contributing to new cell contents now
                   newVolV(piece) = prevMassV(piece) * inv_new_dens
!?!                   newVolV(piece) = prevEcompV(0,piece) / (rho_o*solnData(EINT_VAR,i,j,k)+einternalAdvectedV)
                   if (newVolV(piece) .LE. 0.0) then
                      kappa(piece) = 1.0
                      prevEcompV(0:3,piece) = 0.0
                      newEcompV(0:3,piece) = 0.0
                   else
                      kappa(piece) = prevVolV(piece) / newVolV(piece) !compression factor
                      newEcompV(0:3,piece) = prevEcompV(0:3,piece) * kappa(piece)**(gamComp(0:3) - 1.0)
                      newEinternalV(0:3) = newEinternalV(0:3) + newEcompV(0:3,piece)
                   end if
                end do
             end if

             if(hy_3Ttry_B==1) then 
              aux1 = -dtdx(i) * (oneflx2 - oneflx1)
              if(hy_3Ttry_G==1) then
                 eEleAdvectedPlusPdVV = eEleAdvectedV + 0.5*(p2av1+p2av2)*aux1
                 eIonAdvectedPlusPdVV = eIonAdvectedV + 0.5*(p1av1+p1av2)*aux1
                 eRadAdvectedPlusPdVV = eRadAdvectedV + 0.5*(p3av1+p3av2)*aux1
              else
                 eEleAdvectedPlusPdVV = eEleAdvectedV + solnData(PELE_VAR,i,j,k)*aux1
                 eIonAdvectedPlusPdVV = eIonAdvectedV + solnData(PION_VAR,i,j,k)*aux1
                 eRadAdvectedPlusPdVV = eRadAdvectedV + solnData(PRAD_VAR,i,j,k)*aux1
              end if
             else if(hy_3Ttry_B_rad==1) then
              aux1 = -dtdx(i) * (oneflx2 - oneflx1)
              if(hy_3Ttry_G==1) then
                 eRadAdvectedPlusPdVV = eRadAdvectedV + 0.5*(p3av1+p3av2)*aux1
              else
                 eRadAdvectedPlusPdVV = eRadAdvectedV + solnData(PRAD_VAR,i,j,k)*aux1
              end if
             end if
             if(hy_3Ttry_B==2) then 
              aux1 = -dtdx(i) * ( (eint1flx2 - eint1flx1)               &
                     -0.5e0*(velx_o + solnData(VELX_VAR,i,j,k))        &
                     *aold_t * (p1av2 - p1av1) )
              eIonAdvectedPlusPdVV =  aux1
              aux1 = -dtdx(i) * ( (eint2flx2 - eint2flx1)               &
                     -0.5e0*(velx_o + solnData(VELX_VAR,i,j,k))        &
                     *aold_t * (p2av2 - p2av1) )
              eEleAdvectedPlusPdVV =  aux1
             end if
             if(hy_3Ttry_B_rad==2) then
              aux1 = -dtdx(i) * ( (eint3flx2 - eint3flx1)               &
                     -0.5e0*(velx_o + solnData(VELX_VAR,i,j,k))        &
                     *aold_t * (p3av2 - p3av1) )
              eRadAdvectedPlusPdVV =  aux1
             end if
#define yuckFactor (1.0)
             if(hy_3Ttry_B_rad==3) then
                eRadAdvectedPlusPdVV =  newEinternalV(3) - rho_o*solnData(ERAD_VAR,i,j,k)
             end if
             if(hy_3Ttry_B_rad==0) then
              rescaledERadAdvectedPlusPdVV = eRadAdvectedV
             else
              rescaledERadAdvectedPlusPdVV = eRadAdvectedPlusPdVV
             endif
             if(hy_3Ttry_B==3) then 
                eIonAdvectedPlusPdVV =  newEinternalV(1) - rho_o*solnData(EION_VAR,i,j,k)
                eEleAdvectedPlusPdVV =  newEinternalV(2) - rho_o*solnData(EELE_VAR,i,j,k)
                if ( eintIncreaseV - (eEleAdvectedPlusPdVV+rescaledERadAdvectedPlusPdVV)  &
                     < (eIonAdvectedV) + &
                         (eintIncreaseV - (eIonAdvectedV+eEleAdvectedV+eRadAdvectedV)) &
                           * (solnData(PION_VAR,i,j,k) / solnData(PRES_VAR,i,j,k)) &
                     ) then
#ifdef DEBUG_APR2012
                   print*,' >  ',i,' S:',(.not. forceRageLike),' -> f', &
                           eEleAdvectedPlusPdVV+rescaledERadAdvectedPlusPdVV,   &
                     '>?', eintIncreaseV, solnData(DENS_VAR,i,j,k)*einternal, &
                     solnData(DENS_VAR,i,j,k)/rho_o!,kappa(1),kappa(0),kappa(2)
#endif
                   forceRageLike = .TRUE.
                else if (eIonAdvectedPlusPdVV*yuckFactor  &
                        + eEleAdvectedPlusPdVV+rescaledERadAdvectedPlusPdVV   &
                     > eintIncreaseV) then
#ifdef DEBUG_APR2012
                   print*,' >  ',i,' S:',(.not. forceRageLike),' -> F', &
                        eIonAdvectedPlusPdVV*yuckFactor  &
                        + eEleAdvectedPlusPdVV+rescaledERadAdvectedPlusPdVV,   &
                     ' >', eintIncreaseV, solnData(DENS_VAR,i,j,k)*einternal, &
                     solnData(DENS_VAR,i,j,k)/rho_o!,kappa(1),kappa(0),kappa(2)
#endif
                   forceRageLike = .TRUE.
                else if ( eEleAdvectedPlusPdVV+rescaledERadAdvectedPlusPdVV   &
                     > eintIncreaseV) then
#ifdef DEBUG_APR2012
                   print*,' >  ',i,' S:',(.not. forceRageLike),'0-> F', &
                          eEleAdvectedPlusPdVV+rescaledERadAdvectedPlusPdVV,   &
                     ' >', eintIncreaseV, solnData(DENS_VAR,i,j,k)*einternal, &
                     solnData(DENS_VAR,i,j,k)/rho_o!,kappa(1),kappa(0),kappa(2)
#endif
                   forceRageLike = .TRUE.
                else
#ifdef DEBUG_APR2012
                   print*,'.LE.',i,' S:',(.not. forceRageLike),' keep', &
                          eEleAdvectedPlusPdVV+rescaledERadAdvectedPlusPdVV,   &
                     '<=', eintIncreaseV, solnData(DENS_VAR,i,j,k)*einternal, &
                     solnData(DENS_VAR,i,j,k)/rho_o!,kappa(1),kappa(0),kappa(2)
#endif
                end if
             end if
!!              print*,'eEle:',eEleAdvectedV,eEleAdvectedPlusPdVV
           if((.not. forceRageLike) .AND.(hy_3Ttry_E==2 .OR. &
                (hy_3Ttry_D/=2.0 .AND. hy_3Ttry_D/=3.0))) then !recalibrate based on energy ratio
             if(hy_3Ttry_B==0) then
              rescaledEEleAdvectedPlusPdVV = eEleAdvectedV
              rescaledEIonAdvectedPlusPdVV = eIonAdvectedV
             else
              rescaledEEleAdvectedPlusPdVV = eEleAdvectedPlusPdVV
              rescaledEIonAdvectedPlusPdVV = eIonAdvectedPlusPdVV
             endif
              zeroDummy = 0.0
            if(hy_3Ttry_D/=0.0) then
             if(hy_3Ttry_D==1.0) then
              call Hydro_recalibrateEints(eintIncreaseV-rescaledEEleAdvectedPlusPdVV-rescaledERadAdvectedPlusPdVV,&
                   rescaledEIonAdvectedPlusPdVV,zeroDummy)
             end if
             if(hy_3Ttry_D==1.5) then
              call Hydro_recalibrateEints(eintIncreaseV-rescaledEIonAdvectedPlusPdVV-rescaledERadAdvectedPlusPdVV,&
                   zeroDummy,rescaledEEleAdvectedPlusPdVV)
             end if
             if(hy_3Ttry_D==1.25) then
                if (eintIncreaseV .GE. &
                       rescaledEEleAdvectedPlusPdVV+rescaledERadAdvectedPlusPdVV) then
                   call Hydro_recalibrateEints(eintIncreaseV-rescaledEEleAdvectedPlusPdVV-rescaledERadAdvectedPlusPdVV,&
                        rescaledEIonAdvectedPlusPdVV,zeroDummy)
                else
                   call Hydro_recalibrateEints(eintIncreaseV,&
                        rescaledEIonAdvectedPlusPdVV,rescaledEEleAdvectedPlusPdVV,rescaledERadAdvectedPlusPdVV)
                end if
             end if
             if(hy_3Ttry_D==1.75) then
              call Hydro_recalibrateEints(eintIncreaseV-rescaledEEleAdvectedPlusPdVV,&
                   rescaledEIonAdvectedPlusPdVV,zeroDummy,rescaledERadAdvectedPlusPdVV)
             end if
             if(hy_3Ttry_D==1.875) then
              call Hydro_recalibrateEints(eintIncreaseV-rescaledEIonAdvectedPlusPdVV,&
                   zeroDummy,rescaledEEleAdvectedPlusPdVV,rescaledERadAdvectedPlusPdVV)
             end if
             if(hy_3Ttry_D==2.0) then
             if(hy_3Ttry_F==2) then
              call Hydro_recalibrateEints(eintIncreaseV,&
                   rescaledEIonAdvectedPlusPdVV,rescaledEEleAdvectedPlusPdVV,rescaledERadAdvectedPlusPdVV)
             end if
             end if
             if(hy_3Ttry_D==3.0) then
             if(hy_3Ttry_F==2) then
              call Hydro_recalibrateEints(eintIncreaseV-rescaledERadAdvectedPlusPdVV,&
                   rescaledEIonAdvectedPlusPdVV,rescaledEEleAdvectedPlusPdVV,zeroDummy)
             end if
             end if
            end if
           else if(hy_3Ttry_E==1 .OR. forceRageLike) then  !recalibrate based on pressure ratio
              zeroDummy = 0.0
              PeP = solnData(PELE_VAR,i,j,k) / solnData(PRES_VAR,i,j,k)
              PiP = solnData(PION_VAR,i,j,k) / solnData(PRES_VAR,i,j,k)
              if (hy_3Ttry_D==3.0) then
                 PrP = 0.0
              else
                 PrP = solnData(PRAD_VAR,i,j,k) / solnData(PRES_VAR,i,j,k)
              end if
              eintIncreaseAboveAdvectedV = eintIncreaseV - (eIonAdvectedV+eEleAdvectedV+eRadAdvectedV)
              call Hydro_recalibrateEints(eintIncreaseAboveAdvectedV,&
                   PiP,PeP,PrP)
              rescaledEEleAdvectedPlusPdVV = eEleAdvectedV + PeP
              rescaledEIonAdvectedPlusPdVV = eIonAdvectedV + PiP
              if (hy_3Ttry_D .NE. 3.0 .OR. hy_3Ttry_B_rad==0) then
                 rescaledERadAdvectedPlusPdVV = eRadAdvectedV + PrP
              else
                 rescaledERadAdvectedPlusPdVV = eRadAdvectedPlusPdVV
              end if
              call internal_shiftEints(rescaledEIonAdvectedPlusPdVV, &
                                    rescaledEEleAdvectedPlusPdVV, &
                                    rescaledERadAdvectedPlusPdVV, &
                                    solnData(EION_VAR,i,j,k),solnData(EELE_VAR,i,j,k), &
                                    solnData(ERAD_VAR,i,j,k), &
                                    rho_o, solnData(DENS_VAR,i,j,k),inv_new_dens,i,j,k)

           end if
#ifdef XXHYDRO
              print 800,'eEle:',eEleAdvectedV,eEleAdvectedPlusPdVV,rescaledEEleAdvectedPlusPdVV,&
                   eIonAdvectedV,eIonAdvectedPlusPdVV,rescaledEIonAdvectedPlusPdVV
#endif

!!$              aux1 = -dtdx(i) * ( (eint1flx2 - eint1flx1)               &
!!$                     -0.5e0*(velx_o + solnData(VELX_VAR,i,j,k))        &
!!$                     *aold_t * (p1av2 - p1av1) )
!!$              eIonAdvectedPlusPdVV =  aux1
!!$              aux1 = -dtdx(i) * ( (eint2flx2 - eint2flx1)               &
!!$                     -0.5e0*(velx_o + solnData(VELX_VAR,i,j,k))        &
!!$                     *aold_t * (p2av2 - p2av1) )
!!$              eEleAdvectedPlusPdVV =  aux1
!!$              aux1 = -dtdx(i) * ( (eint3flx2 - eint3flx1)               &
!!$                     -0.5e0*(velx_o + solnData(VELX_VAR,i,j,k))        &
!!$                     *aold_t * (p3av2 - p3av1) )
!!$              eRadAdvectedPlusPdVV =  aux1
!!$              rescaledEEleAdvectedPlusPdVV = eEleAdvectedPlusPdVV
!!$              rescaledEIonAdvectedPlusPdVV = eIonAdvectedPlusPdVV
!!$              rescaledERadAdvectedPlusPdVV = eRadAdvectedPlusPdVV
!!$              zeroDummy = 0.0
!!$              call Hydro_recalibrateEints(eintIncreaseV-rescaledEEleAdvectedPlusPdVV,&
!!$                   rescaledEIonAdvectedPlusPdVV,zeroDummy,rescaledERadAdvectedPlusPdVV)
!!$!!              call Hydro_recalibrateEints(eintIncreaseV,&
!!$!!                   rescaledEIonAdvectedPlusPdVV,rescaledEEleAdvectedPlusPdVV,rescaledERadAdvectedPlusPdVV)
!!$              print 800,'eEle+',eEleAdvectedV,eEleAdvectedPlusPdVV,rescaledEEleAdvectedPlusPdVV,&
!!$                   eIonAdvectedV,eIonAdvectedPlusPdVV,rescaledEIonAdvectedPlusPdVV


!********** E1_VAR BEGIN **********
#ifdef FLASH_MULTISPECIES
              call Eos_getAbarZbar(solnVec=solnData(:,i,j,k),Ye=Ye)           !DEV: for now
#else
#if defined(SUMY_MSCALAR) && defined(YE_MSCALAR)
              Ye = solnData(YE_MSCALAR,i,j,k) !DEV: should not matter - to be removed? - KW
#else
              call Eos_getAbarZbar(solnVec=solnData(:,i,j,k),Ye=Ye)           !DEV: for now
#endif
#endif
              ekinElecFrac = Ye * hy_eMassInUAmu
              ekinIonFrac = 1 - ekinElecFrac
              ! update the ion energy
              aux1 = - dtdx(i) * (e1flx2 - e1flx1)                             &
                   + dt*0.5e00 * (1 - eKinElecFrac)                            &
                   *( rho_o*velx_o*grav1d_o                                  &
                     +solnData(DENS_VAR,i,j,k)*solnData(VELX_VAR,i,j,k)*grav1d &
                    )

#ifdef E1_VAR
              etot = (rho_o * solnData(  E1_VAR,i,j,k) + aux1)*inv_new_dens
#endif
              
#ifdef EION_VAR
              ! get the ion internal energy
              einternal = solnData(EION_VAR,i,j,k)

              ! update ion internal energy
              aux1 = -dtdx(i) * ( (eint1flx2 - eint1flx1)               &
                     -0.5e0*(velx_o + solnData(VELX_VAR,i,j,k))        &
                     *aold_t * (p1av2 - p1av1) )
              
#ifdef DEBUG_MAR2012
              if(rho_o * einternal + rescaledEIonAdvectedPlusPdVV .LE.0.0) then
!!$                 print*,'NegION! at',i,(rho_o * einternal),rescaledEIonAdvectedPlusPdVV
                 print*,'NegION! at',i,einternal,(rho_o * einternal*inv_new_dens),rescaledEIonAdvectedPlusPdVV*inv_new_dens
              elseif( rescaledEIonAdvectedPlusPdVV / (rho_o * einternal) > 10.0) then
!!$                 print*,'PosION: at',i,(rho_o * einternal),rescaledEIonAdvectedPlusPdVV
                 print*,'PosION: at',i,einternal,(rho_o * einternal*inv_new_dens),rescaledEIonAdvectedPlusPdVV*inv_new_dens
              end if
#endif
!              einternal = (rho_o * einternal + aux1)*inv_new_dens
              einternal = (rho_o * einternal + rescaledEIonAdvectedPlusPdVV)*inv_new_dens
!!$              einternal = max(einternal, hy_smallp*inv_new_dens)

! compute the new kinetic energy       
              ekin = 0.5e0 * ( solnData(VELX_VAR,i,j,k)**2             &
                              +solnData(VELY_VAR,i,j,k)**2             &
                              +solnData(VELZ_VAR,i,j,k)**2)
              ekin = ekin * (1 - ekinElecFrac)

! test whether we should use the internal energy from the evolution
              if (.TRUE. .OR. einternal .LT. hy_eint1Switch*ekin) then
                 etot = einternal + ekin
              else
                 einternal = etot - ekin
              endif
              
              if (inShock) then
                 eintIncreaseV = solnData(DENS_VAR,i,j,k)*einternal - rho_o*solnData(EION_VAR,i,j,k)
                 ! update internal energy advected
                 aux1 = -dtdx(i) * (eia1flx2 - eia1flx1)

                 einternalAdvectedV = aux1

                 eionIncreaseAboveAdvectedV =  eintIncreaseV - einternalAdvectedV
#ifdef DEBUG_XHYDRO
998              format(1x,I2,' inShock Ions:',3(1x,1PG22.15))
                 print 998,i,eintIncreaseV,einternalAdvectedV,eionIncreaseAboveAdvectedV
#endif
!!$                 eionTargetIncreaseV = max(eionIncreaseAboveAdvectedV, ekinIonFrac*eintIncreaseAboveAdvectedV)
!!$                 eeleTargetIncreaseV = max(eionIncreaseAboveAdvectedV, ekinEleFrac*eintIncreaseAboveAdvectedV)
!!$                 eionAdjustment = max(0.0,eionTargetIncreaseV-eionIncreaseAboveAdvectedV) * inv_new_dens
!!$                 if (eionAdjustment .LE. 0.0) then
!!$                    print*,i,' InShock ignoring negative eionAdjustment:',eionAdjustment
!!$                    eionAdjustment = 0.0
!!$                    inShock = .FALSE.
!!$                 else
!!$                    einternal = einternal + eionAdjustment
!!$                    etot = etot + eionAdjustment
!!$                 end if
              end if

#ifdef E1_VAR
              solnData(E1_VAR,i,j,k) = etot
#endif
              solnData(EION_VAR,i,j,k) = einternal
#else
              
#ifdef E1_VAR
              solnData(E1_VAR,i,j,k) = etot
#endif
#endif
!********** E1_VAR END **********

!********** E2_VAR BEGIN **********
#ifdef FLASH_MULTISPECIES
              call Eos_getAbarZbar(solnVec=solnData(:,i,j,k),Ye=Ye)           !DEV: for now
#else
#if defined(SUMY_MSCALAR) && defined(YE_MSCALAR)
              Ye = solnData(YE_MSCALAR,i,j,k) !DEV: should not matter - to be removed? - KW
#else
              call Eos_getAbarZbar(solnVec=solnData(:,i,j,k),Ye=Ye)           !DEV: for now
#endif
#endif
              ekinElecFrac = Ye * hy_eMassInUAmu
              ! update the electron energy
              aux1 = - dtdx(i) * (e2flx2 - e2flx1)                             &
                   + dt*0.5e00 * ekinElecFrac                                &
                   *( rho_o*velx_o*grav1d_o                                  &
                     +solnData(DENS_VAR,i,j,k)*solnData(VELX_VAR,i,j,k)*grav1d &
                    )

#ifdef E2_VAR
              etot = (rho_o * solnData(  E2_VAR,i,j,k) + aux1)*inv_new_dens
#endif
              
#ifdef EELE_VAR
              ! get the electron internal energy
              einternal = solnData(EELE_VAR,i,j,k)

              ! update electron internal energy
              aux1 = -dtdx(i) * ( (eint2flx2 - eint2flx1)               &
                     -0.5e0*(velx_o + solnData(VELX_VAR,i,j,k))        &
                     *aold_t * (p2av2 - p2av1) )
              
!              einternal = (rho_o * einternal + aux1)*inv_new_dens
              einternal = (rho_o * einternal + rescaledEEleAdvectedPlusPdVV)*inv_new_dens
!!$              einternal = max(einternal, hy_smallp*inv_new_dens)

! compute the new kinetic energy       
              ekin = 0.5e0 * ( solnData(VELX_VAR,i,j,k)**2             &
                              +solnData(VELY_VAR,i,j,k)**2             &
                              +solnData(VELZ_VAR,i,j,k)**2)
              ekin = ekin * ekinElecFrac

! test whether we should use the internal energy from the evolution
              if (.TRUE. .OR. einternal .LT. hy_eint2Switch*ekin) then
                 etot = einternal + ekin
              else
                 einternal = etot - ekin
              endif
              
              if (inShock .AND. .FALSE.) then
                 eionTargetIncreaseV = max(eionIncreaseAboveAdvectedV, ekinIonFrac*eintIncreaseAboveAdvectedV)
                 eionTargetIncreaseV = min(eionTargetIncreaseV, einternal - hy_smallp)
!!$                 eeleTargetIncreaseV = max(eionIncreaseAboveAdvectedV, ekinEleFrac*eintIncreaseAboveAdvectedV)
                 eionAdjustment = max(0.0,eionTargetIncreaseV-eionIncreaseAboveAdvectedV) * inv_new_dens
                 if (eionAdjustment .LE. 0.0) then
#ifdef DEBUG_XHYDRO
                    print*,i,' InShock ignoring negative eionAdjustment:',eionAdjustment
#endif
                    eionAdjustment = 0.0
                    inShock = .FALSE.
                 else
                    solnData(EION_VAR,i,j,k) = solnData(EION_VAR,i,j,k) + eionAdjustment
#ifdef E1_VAR
                    solnData(E1_VAR,i,j,k) = solnData(E1_VAR,i,j,k) + eionAdjustment
#endif
                    eeleAdjustment = - min(eionAdjustment, einternal)
#ifdef DEBUG_XHYDRO
                    if (eionAdjustment + eeleAdjustment > 0.0) then
                       print*,i,' InShock WARNING some of adjustments:',eionAdjustment+eeleAdjustment
                    end if
#endif
                    einternal = einternal + eeleAdjustment
                    if (einternal < hy_smallp*inv_new_dens) then
                       print*,i,' InShock WARNING electron energy now very small:',einternal
                    end if
                    etot = etot + eeleAdjustment
                 end if
              end if

#ifdef E2_VAR
              solnData(E2_VAR,i,j,k) = etot
#endif
              solnData(EELE_VAR,i,j,k) = einternal
#else
              
#ifdef E2_VAR
              solnData(E2_VAR,i,j,k) = etot
#endif
#endif
!********** E2_VAR END **********

!********** E3_VAR BEGIN **********
              ! update the radiation energy
              aux1 = - dtdx(i) * (e3flx2 - e3flx1)                             &
                   + dt*0.5e00  * 0.0                                              &
                   *( rho_o*velx_o*grav1d_o                                  &
                     +solnData(DENS_VAR,i,j,k)*solnData(VELX_VAR,i,j,k)*grav1d &
                    )

#ifdef E3_VAR
              etot = (rho_o * solnData(  E3_VAR,i,j,k) + aux1)*inv_new_dens
#endif
              
#ifdef ERAD_VAR
              ! get the radiation internal energy
              einternal = solnData(ERAD_VAR,i,j,k)

              ! update radiation internal energy
              aux1 = -dtdx(i) * ( (eint3flx2 - eint3flx1)               &
                     -0.5e0*(velx_o + solnData(VELX_VAR,i,j,k))        &
                     *aold_t * (p3av2 - p3av1) )
              
!              einternal = (rho_o * einternal + aux1)*inv_new_dens
              einternal = (rho_o * einternal + rescaledERadAdvectedPlusPdVV)*inv_new_dens
!!$              einternal = max(einternal, hy_smallp*inv_new_dens) !! ERAD can be zero in 2T mode.

! compute the new kinetic energy        - ZERO for now (and perhaps forever)
              ekin = 0.0

! test whether we should use the internal energy from the evolution
!!$              if (einternal .LT. hy_eint3Switch*ekin) then
                 etot = einternal + ekin
!!$              else
!!$                 einternal = etot - ekin
!!$              endif
              
#ifdef E3_VAR
              solnData(E3_VAR,i,j,k) = etot
#endif
              solnData(ERAD_VAR,i,j,k) = einternal
#else
              
#ifdef E3_VAR
              solnData(E3_VAR,i,j,k) = etot
#endif
#endif
!********** E3_VAR END **********

              call internal_shiftEints(solnData(EION_VAR,i,j,k), &
                                    solnData(EELE_VAR,i,j,k), &
                                    solnData(ERAD_VAR,i,j,k), &
                                    0.0,0.0,0.0, &
                                    1.,1.,1.,i,j,k)
           end do
        end do
     end do
     
  end if
               
!===============================================================================

! update in y direction
  
  if ((NDIM >= 2) .and. (xyzswp == SWEEP_Y)) then

     select case (rangeSwitch)
     case (UPDATE_INTERIOR)
        jmin  = blkLimits(LOW,JAXIS) + 1 ! NGUARD*K2D + 2
        jmax  = blkLimits(HIGH,JAXIS) - 1 ! NGUARD*K2D + NYB - 1
        jskip = 1
     case (UPDATE_BOUND)
        jmin  = blkLimits(LOW,JAXIS) ! NGUARD*K2D + 1
        jmax  = blkLimits(HIGH,JAXIS) ! NGUARD*K2D + NYB
        jskip = blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) !NYB - 1
     case default
        jmin  = blkLimits(LOW,JAXIS) ! NGUARD*K2D + 1
        jmax  = blkLimits(HIGH,JAXIS) ! NGUARD*K2D + NYB
        jskip = 1
     end select

#ifdef DEBUG_HYDR
     print*,'the sweep direction is',xyzswp      ! within DEBUG
     print*,'the blkLimits is',blkLimits         ! within DEBUG
     print*,'the blkLimitsGC is',blkLimitsGC     ! within DEBUG
     print*,'the update mode is',rangeSwitch,' numCells',numCells   ! within DEBUG
     print*,'jmin etc',jmin,jmax,jskip           ! within DEBUG
#endif

     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = jmin, jmax, jskip
           do i = blkLimits(LOW,IAXIS) ,blkLimits(HIGH,IAXIS)

              dtdx(j) = tempDtDx(i,j,k)

              rhoflx1  = tempFlx(RHO_FLUX,i,j,k)
              rhoflx2  = tempFlx(RHO_FLUX,i,j+1,k)

              uflx1    = tempFlx(U_FLUX,i,j,k)
              uflx2    = tempFlx(U_FLUX,i,j+1,k)

              pav1     = tempFlx(P_FLUX,i,j,k)
              pav2     = tempFlx(P_FLUX,i,j+1,k)

              p1av1     = tempFlx(PION_FLUX,i,j,k)
              p1av2     = tempFlx(PION_FLUX,i,j+1,k)
              p2av1     = tempFlx(PELE_FLUX,i,j,k)
              p2av2     = tempFlx(PELE_FLUX,i,j+1,k)
              p3av1     = tempFlx(PRAD_FLUX,i,j,k)
              p3av2     = tempFlx(PRAD_FLUX,i,j+1,k)

              utflx1   = tempFlx(UT_FLUX,i,j,k)
              utflx2   = tempFlx(UT_FLUX,i,j+1,k)

              uttflx1  = tempFlx(UTT_FLUX,i,j,k)
              uttflx2  = tempFlx(UTT_FLUX,i,j+1,k)

              eflx1    = tempFlx(E_FLUX,i,j,k)
              eflx2    = tempFlx(E_FLUX,i,j+1,k)

              e1flx1    = tempFlx(E1_FLUX,i,j,k)
              e1flx2    = tempFlx(E1_FLUX,i,j+1,k)
              e2flx1    = tempFlx(E2_FLUX,i,j,k)
              e2flx2    = tempFlx(E2_FLUX,i,j+1,k)
              e3flx1    = tempFlx(E3_FLUX,i,j,k)
              e3flx2    = tempFlx(E3_FLUX,i,j+1,k)

              eintflx1 = tempFlx(EINT_FLUX,i,j,k)
              eintflx2 = tempFlx(EINT_FLUX,i,j+1,k)

              eint1flx1 = tempFlx(EION_FLUX,i,j,k)
              eint1flx2 = tempFlx(EION_FLUX,i,j+1,k)
              eint2flx1 = tempFlx(EELE_FLUX,i,j,k)
              eint2flx2 = tempFlx(EELE_FLUX,i,j+1,k)
              eint3flx1 = tempFlx(ERAD_FLUX,i,j,k)
              eint3flx2 = tempFlx(ERAD_FLUX,i,j+1,k)

              eiaflx1 = tempFlx(EIA_FLUX,i,j,k)
              eiaflx2 = tempFlx(EIA_FLUX,i,j+1,k)

              eia1flx1 = tempFlx(EI1A_FLUX,i,j,k)
              eia1flx2 = tempFlx(EI1A_FLUX,i,j+1,k)
              eia2flx1 = tempFlx(EI2A_FLUX,i,j,k)
              eia2flx2 = tempFlx(EI2A_FLUX,i,j+1,k)
              eia3flx1 = tempFlx(EI3A_FLUX,i,j,k)
              eia3flx2 = tempFlx(EI3A_FLUX,i,j+1,k)

              oneflx1 = tempFlx(ONE_FLUX,i,j,k)
              oneflx2 = tempFlx(ONE_FLUX,i,j+1,k)
#ifdef VOLD_FLUX
              voFlx(1) =  tempFlx(VOLD_FLUX,i,j,k)
              voFlx(2) = -tempFlx(VOLD_FLUX,i,j+1,k)
#endif

              do kk = 1,hy_numXn
                 xnflx1(kk) = tempFlx(SPECIES_FLUX_BEGIN+kk-1,i,j,k)
                 xnflx2(kk) = tempFlx(SPECIES_FLUX_BEGIN+kk-1,i,j+1,k)
              end do

              aold_t   = tempArea(i,j,k)
              grav1d_o = tempGrav1d_o(i,j,k)
              grav1d   = tempGrav1d(i,j,k)
              fict1d   = tempFict(i,j,k)

              rho_o    = solnData(DENS_VAR,i,j,k)
              prevMassV(0) = rho_o
              vely_o   = solnData(VELY_VAR,i,j,k)

              if ( hy_useCmaAdvection .and. NSPECIES > 1  ) then

! update the partial mass densities and passive scalars * density

                 prevMassV(1) = 0.0; prevMassV(2) = 0.0
                 do n = 1, hy_numXn
                    solnData(SPECIES_BEGIN-1+n,i,j,k) =               &
                              rho_o*solnData(SPECIES_BEGIN-1+n,i,j,k) &
                            -dtdx(j)*(xnflx2(n) - xnflx1(n))
                    prevMassV(1) = prevMassV(1) + dtdx(j) * xnflx1(n)
                    prevMassV(2) = prevMassV(2) - dtdx(j) * xnflx2(n)
                 end do

! update the total mass density

                 solnData(DENS_VAR,i,j,k) = 0.e0

                 do n = 1, NSPECIES
                    solnData(DENS_VAR,i,j,k) =                        &
                    solnData(DENS_VAR,i,j,k) +                        &
                    solnData(SPECIES_BEGIN-1+n,i,j,k)
                 end do

! recover partial densities and passive scalars * density

                 inv_new_dens = 1.e0/solnData(DENS_VAR,i,j,k)

                 do n = 1, hy_numXn
                    solnData(SPECIES_BEGIN-1+n,i,j,k) =               &
                    solnData(SPECIES_BEGIN-1+n,i,j,k) * inv_new_dens
                 end do

              else

! update the density               
                 solnData(DENS_VAR,i,j,k) = solnData(DENS_VAR,i,j,k) - &
                                       dtdx(j) * (rhoflx2-rhoflx1)
                 prevMassV(1) =  dtdx(j) * rhoflx1
                 prevMassV(2) = -dtdx(j) * rhoflx2

! update the mass fractions and passive scalars
                 do n = 1, hy_numXn
                    solnData(SPECIES_BEGIN-1+n,i,j,k) =               &
                           ( rho_o*solnData(SPECIES_BEGIN-1+n,i,j,k)  &
                            -dtdx(j)*(xnflx2(n) - xnflx1(n))         &
                           )/solnData(DENS_VAR,i,j,k)
                 end do

              end if

! limit the density 
              solnData(DENS_VAR,i,j,k) =                              &
                   max(hy_smlrho, solnData(DENS_VAR,i,j,k))

              inv_new_dens = 1.e0/solnData(DENS_VAR,i,j,k)
!#ifdef CIP
!
!              if ( itrcr > 0 ) then
!
!                 ! CIP update and limit tracer
!
!                 do n = itrcr,itrcr
!                    solnData(n,i,j,k) = atan(solnData(n,i,j,k))/(0.9999d0*pi) + 0.5e0
!                    solnData(n,i,j,k) = max(0.e0, min(1.e0, solnData(n,i,j,k) ))
!                 end do
!
!              end if
!#endif
              
! update the velocities               
              aux1 = -dtdx(j)*( (uflx2 - uflx1)                      &
                               +aold_t*(pav2 - pav1))                &
               +0.5e0*dt                                             &
               *( (solnData(DENS_VAR,i,j,k) + rho_o)*fict1d           &
                 +(rho_o*grav1d_o + solnData(DENS_VAR,i,j,k)*grav1d)  &
                )

              solnData(VELY_VAR,i,j,k) =                              &
                         ( rho_o * solnData(VELY_VAR,i,j,k)           &
                          +aux1                                      &
                         )*inv_new_dens
              
              solnData(VELX_VAR,i,j,k) =                              &
                         ( rho_o * solnData(VELX_VAR,i,j,k)           &
                          -dtdx(j) * (utflx2 - utflx1)               &
                         )*inv_new_dens
              
              solnData(VELZ_VAR,i,j,k) =                              &
                         ( rho_o * solnData(VELZ_VAR,i,j,k)           &
                          -dtdx(j) * (uttflx2 - uttflx1)             &
                         )*inv_new_dens

! update the total energy
              aux1 = - dtdx(j) * (eflx2 - eflx1)                             &
                   + dt*0.5e00                                               &
                   *( rho_o*vely_o*grav1d_o                                  &
                     +solnData(DENS_VAR,i,j,k)*solnData(VELY_VAR,i,j,k)*grav1d &
                    )

              etot = (rho_o * solnData(ENER_VAR,i,j,k) + aux1)*inv_new_dens

#ifdef EINT_VAR
! get the internal energy
              einternal = solnData(EINT_VAR,i,j,k)
              prevEcompV(0,0) = rho_o * solnData(EINT_VAR,i,j,k)
               
! update internal energy
              aux1 = -dtdx(j) * ( (eintflx2 - eintflx1)               &
                     -0.5e0*(vely_o + solnData(VELY_VAR,i,j,k))        &
                     *aold_t * (pav2 - pav1) )

              einternal = (rho_o * einternal + aux1)*inv_new_dens
              einternal = max(einternal, hy_smallp*inv_new_dens)

! compute the new kinetic energy
              ekin = 0.5e0 * ( solnData(VELX_VAR,i,j,k)**2             &
                              +solnData(VELY_VAR,i,j,k)**2             &
                              +solnData(VELZ_VAR,i,j,k)**2)

! test whether we should use the internal energy from the evolution
              if (einternal .LT. hy_eintSwitch*ekin) then
                 etot = einternal + ekin
              else
                 einternal = etot - ekin
              endif
              
#ifdef DBGS_VAR
              solndata(DBGS_VAR,i,j,k) = 0.5 * &
                   ( solndata(DBGS_VAR,i,j,k)  + tempFlx(SHOK_FLUX,i,j,k)+tempFlx(SHOK_FLUX,i,j+1,k) )
#endif
              inShock = ( tempFlx(SHOK_FLUX,i,j,k)+tempFlx(SHOK_FLUX,i,j+1,k) > 0 )
              if (hy_3Ttry_useShockDetect) then
                 forceRageLike = (.NOT. inShock)
              else
                 forceRageLike = .FALSE.
              end if

              eintIncreaseV = solnData(DENS_VAR,i,j,k)*einternal - rho_o*solnData(EINT_VAR,i,j,k)
              ! update internal energy advected
              aux1 = -dtdx(j) * (eiaflx2 - eiaflx1)
              einternalAdvectedV = aux1
              eintIncreaseAboveAdvectedV =  eintIncreaseV - einternalAdvectedV

              solnData(ENER_VAR,i,j,k) = etot
              solnData(EINT_VAR,i,j,k) = einternal
#else
              
              solnData(ENER_VAR,i,j,k) = etot
#endif

!!$              if(hy_3Ttry_B==0 .OR. hy_3Ttry_B==1 .OR. (hy_3Ttry_D .GE. 2.0 .AND. hy_3Ttry_E==1)) then
                 aux1 = -dtdx(j) * (eia3flx2 - eia3flx1)
                 eRadAdvectedV = aux1
                 aux1 = -dtdx(j) * (eia2flx2 - eia2flx1)
                 eEleAdvectedV = aux1
                 aux1 = -dtdx(j) * (eia1flx2 - eia1flx1)
                 eIonAdvectedV = aux1
!!$              else if(hy_3Ttry_B_rad==0 .OR. hy_3Ttry_B_rad==1) then
!!$                 aux1 = -dtdx(j) * (eia3flx2 - eia3flx1)
!!$                 eRadAdvectedV = aux1
!!$              end if

             if(hy_3Ttry_B==3 .OR. hy_3Ttry_B_rad==3) then
                prevVolV(0) = 1.0
                prevVolV(1) = dtdx(j) * voFlx(1)
                prevVolV(2) = dtdx(j) * voFlx(2)
                ! Note prevEcompV(0,0) was initialized above from solnData(EINT_VAR,i,j,k), before that changed.
                prevEcompV(1,0) = rho_o * solnData(EION_VAR,i,j,k)
                prevEcompV(2,0) = rho_o * solnData(EELE_VAR,i,j,k)
                prevEcompV(3,0) = rho_o * solnData(ERAD_VAR,i,j,k)
                do side = 1,2      !first left then right interface
                   select case (side)
                   case(1)
                      prevEcompV(0,side) = dtdx(j) * eiaflx1
                      prevEcompV(1,side) = dtdx(j) * eia1flx1
                      prevEcompV(2,side) = dtdx(j) * eia2flx1
                      prevEcompV(3,side) = dtdx(j) * eia3flx1
                   case(2)
                      prevEcompV(0,side) = - dtdx(j) * eiaflx2
                      prevEcompV(1,side) = - dtdx(j) * eia1flx2
                      prevEcompV(2,side) = - dtdx(j) * eia2flx2
                      prevEcompV(3,side) = - dtdx(j) * eia3flx2
                   end select
                   if (voFlx(side) .GE. 0.0) then ! Stuff is coming in on this side
                      ! ... keep it. 
                   else                           ! Stuff is going out on this side
                      prevMassV(0) = prevMassV(0) + prevMassV(side)
                      prevVolV(0)  = prevVolV(0)  + prevVolV(side)
                      prevEcompV(0:3,0)  = prevEcompV(0:3,0)  + prevEcompV(0:3,side)
                      prevMassV(side) = 0.0
                      prevVolV(side) = 0.0
                   end if
                end do

                newEinternalV(0:3) = 0.0
                do piece = 0,2      !all fluid pieces contributing to new cell contents now
                   newVolV(piece) = prevMassV(piece) * inv_new_dens
                   if (newVolV(piece) .LE. 0.0) then
                      kappa(piece) = 1.0
                      prevEcompV(0:3,piece) = 0.0
                      newEcompV(0:3,piece) = 0.0
                   else
                      kappa(piece) = prevVolV(piece) / newVolV(piece) !compression factor
                      newEcompV(0:3,piece) = prevEcompV(0:3,piece) * kappa(piece)**(gamComp(0:3) - 1.0)
                      newEinternalV(0:3) = newEinternalV(0:3) + newEcompV(0:3,piece)
                   end if
                end do
              end if
 
              if(hy_3Ttry_B==1) then 
                 aux1 = -dtdx(j) * (oneflx2 - oneflx1)
                 if(hy_3Ttry_G==1) then
                    eEleAdvectedPlusPdVV = eEleAdvectedV + 0.5*(p2av1+p2av2)*aux1
                    eIonAdvectedPlusPdVV = eIonAdvectedV + 0.5*(p1av1+p1av2)*aux1
                    eRadAdvectedPlusPdVV = eRadAdvectedV + 0.5*(p3av1+p3av2)*aux1
                 else
                    eEleAdvectedPlusPdVV = eEleAdvectedV + solnData(PELE_VAR,i,j,k)*aux1
                    eIonAdvectedPlusPdVV = eIonAdvectedV + solnData(PION_VAR,i,j,k)*aux1
                    eRadAdvectedPlusPdVV = eRadAdvectedV + solnData(PRAD_VAR,i,j,k)*aux1
                 end if
              else if(hy_3Ttry_B_rad==1) then 
                 aux1 = -dtdx(j) * (oneflx2 - oneflx1)
                 if(hy_3Ttry_G==1) then
                    eRadAdvectedPlusPdVV = eRadAdvectedV + 0.5*(p3av1+p3av2)*aux1
                 else
                    eRadAdvectedPlusPdVV = eRadAdvectedV + solnData(PRAD_VAR,i,j,k)*aux1
                 end if
              end if
              if(hy_3Ttry_B==2) then 
                 aux1 = -dtdx(j) * ( (eint1flx2 - eint1flx1)               &
                      -0.5e0*(vely_o + solnData(VELY_VAR,i,j,k))        &
                      *aold_t * (p1av2 - p1av1) )
                 eIonAdvectedPlusPdVV =  aux1
                 aux1 = -dtdx(j) * ( (eint2flx2 - eint2flx1)               &
                      -0.5e0*(vely_o + solnData(VELY_VAR,i,j,k))        &
                      *aold_t * (p2av2 - p2av1) )
                 eEleAdvectedPlusPdVV =  aux1
              end if
              if(hy_3Ttry_B_rad==2) then 
                 aux1 = -dtdx(j) * ( (eint3flx2 - eint3flx1)               &
                      -0.5e0*(vely_o + solnData(VELY_VAR,i,j,k))        &
                      *aold_t * (p3av2 - p3av1) )
                 eRadAdvectedPlusPdVV =  aux1
              end if
!!#define yuckFactor (1.0)
             if(hy_3Ttry_B_rad==3) then
                eRadAdvectedPlusPdVV =  newEinternalV(3) - rho_o*solnData(ERAD_VAR,i,j,k)
             end if
             if(hy_3Ttry_B_rad==0) then
              rescaledERadAdvectedPlusPdVV = eRadAdvectedV
             else
              rescaledERadAdvectedPlusPdVV = eRadAdvectedPlusPdVV
             endif
             if(hy_3Ttry_B==3) then 
                eIonAdvectedPlusPdVV =  newEinternalV(1) - rho_o*solnData(EION_VAR,i,j,k)
                eEleAdvectedPlusPdVV =  newEinternalV(2) - rho_o*solnData(EELE_VAR,i,j,k)
                if ( eintIncreaseV - (eEleAdvectedPlusPdVV+rescaledERadAdvectedPlusPdVV)  &
                     < (eIonAdvectedV) + &
                         (eintIncreaseV - (eIonAdvectedV+eEleAdvectedV+eRadAdvectedV)) &
                           * (solnData(PION_VAR,i,j,k) / solnData(PRES_VAR,i,j,k)) &
                     ) then
                   forceRageLike = .TRUE.
                else if (eIonAdvectedPlusPdVV*yuckFactor  &
                        + eEleAdvectedPlusPdVV+rescaledERadAdvectedPlusPdVV   &
                     > eintIncreaseV) then
                   forceRageLike = .TRUE.
                else if ( eEleAdvectedPlusPdVV+rescaledERadAdvectedPlusPdVV   &
                     > eintIncreaseV) then
                   forceRageLike = .TRUE.
                end if
             end if
              !!              print*,'eEleY',eEleAdvectedV,eEleAdvectedPlusPdVV
              if((.not. forceRageLike) .AND.(hy_3Ttry_E==2 .OR. &
                   ((hy_3Ttry_D/=2.0).AND.(hy_3Ttry_D/=3.0)))) then !recalibrate based on energy ratio
                 if(hy_3Ttry_B==0) then
                    rescaledEEleAdvectedPlusPdVV = eEleAdvectedV
                    rescaledEIonAdvectedPlusPdVV = eIonAdvectedV
                 else
                    rescaledEEleAdvectedPlusPdVV = eEleAdvectedPlusPdVV
                    rescaledEIonAdvectedPlusPdVV = eIonAdvectedPlusPdVV
                 endif
                 zeroDummy = 0.0
                 if(hy_3Ttry_D/=0.0) then
                    if(hy_3Ttry_D==1.0) then
                       call Hydro_recalibrateEints(eintIncreaseV-rescaledEEleAdvectedPlusPdVV-rescaledERadAdvectedPlusPdVV,&
                            rescaledEIonAdvectedPlusPdVV,zeroDummy)
                    end if
                    if(hy_3Ttry_D==1.5) then
                       call Hydro_recalibrateEints(eintIncreaseV-rescaledEIonAdvectedPlusPdVV-rescaledERadAdvectedPlusPdVV,&
                            zeroDummy,rescaledEEleAdvectedPlusPdVV)
                    end if
                    if(hy_3Ttry_D==1.25) then
                       if (eintIncreaseV .GE. &
                            rescaledEEleAdvectedPlusPdVV+rescaledERadAdvectedPlusPdVV) then
                          call Hydro_recalibrateEints(eintIncreaseV-rescaledEEleAdvectedPlusPdVV-rescaledERadAdvectedPlusPdVV,&
                               rescaledEIonAdvectedPlusPdVV,zeroDummy)
                       else
                          call Hydro_recalibrateEints(eintIncreaseV,&
                               rescaledEIonAdvectedPlusPdVV,rescaledEEleAdvectedPlusPdVV,rescaledERadAdvectedPlusPdVV)
                       end if
                    end if
                    if(hy_3Ttry_D==1.75) then
                       call Hydro_recalibrateEints(eintIncreaseV-rescaledEEleAdvectedPlusPdVV,&
                            rescaledEIonAdvectedPlusPdVV,zeroDummy,rescaledERadAdvectedPlusPdVV)
                    end if
                    if(hy_3Ttry_D==1.875) then
                       call Hydro_recalibrateEints(eintIncreaseV-rescaledEIonAdvectedPlusPdVV,&
                            zeroDummy,rescaledEEleAdvectedPlusPdVV,rescaledERadAdvectedPlusPdVV)
                    end if
                    if(hy_3Ttry_D==2.0) then
                       if(hy_3Ttry_F==2) then
                          call Hydro_recalibrateEints(eintIncreaseV,&
                               rescaledEIonAdvectedPlusPdVV,rescaledEEleAdvectedPlusPdVV,rescaledERadAdvectedPlusPdVV)
                       end if
                    end if
                    if(hy_3Ttry_D==3.0) then
                       if(hy_3Ttry_F==2) then
                          call Hydro_recalibrateEints(eintIncreaseV-rescaledERadAdvectedPlusPdVV,&
                               rescaledEIonAdvectedPlusPdVV,rescaledEEleAdvectedPlusPdVV,zeroDummy)
                       end if
                    end if
                 end if
              else if(hy_3Ttry_E==1 .OR. forceRageLike) then  !recalibrate based on pressure ratio
                 zeroDummy = 0.0
                 PeP = solnData(PELE_VAR,i,j,k) / solnData(PRES_VAR,i,j,k)
                 PiP = solnData(PION_VAR,i,j,k) / solnData(PRES_VAR,i,j,k)
                 if (hy_3Ttry_D==3.0) then
                    PrP = 0.0
                 else
                    PrP = solnData(PRAD_VAR,i,j,k) / solnData(PRES_VAR,i,j,k)
                 end if
                 eintIncreaseAboveAdvectedV = eintIncreaseV - (eIonAdvectedV+eEleAdvectedV+eRadAdvectedV)
                 call Hydro_recalibrateEints(eintIncreaseAboveAdvectedV,&
                      PiP,PeP,PrP)
                 rescaledEEleAdvectedPlusPdVV = eEleAdvectedV + PeP
                 rescaledEIonAdvectedPlusPdVV = eIonAdvectedV + PiP
                 if (hy_3Ttry_D .NE. 3.0 .OR. hy_3Ttry_B_rad==0) then
                    rescaledERadAdvectedPlusPdVV = eRadAdvectedV + PrP
                 else
                    rescaledERadAdvectedPlusPdVV = eRadAdvectedPlusPdVV
                 end if
                 call internal_shiftEints(rescaledEIonAdvectedPlusPdVV, &
                                       rescaledEEleAdvectedPlusPdVV, &
                                       rescaledERadAdvectedPlusPdVV, &
                                       solnData(EION_VAR,i,j,k),solnData(EELE_VAR,i,j,k), &
                                       solnData(ERAD_VAR,i,j,k), &
                                       rho_o, solnData(DENS_VAR,i,j,k),inv_new_dens,i,j,k)
              end if
#ifdef XXHYDRO
              print 800,'eEleY',eEleAdvectedV,eEleAdvectedPlusPdVV,rescaledEEleAdvectedPlusPdVV,&
                   eIonAdvectedV,eIonAdvectedPlusPdVV,rescaledEIonAdvectedPlusPdVV
#endif


!********** E1_VAR BEGIN **********
#ifdef FLASH_MULTISPECIES
              call Eos_getAbarZbar(solnVec=solnData(:,i,j,k),Ye=Ye)           !DEV: for now
#else
#if defined(SUMY_MSCALAR) && defined(YE_MSCALAR)
              Ye = solnData(YE_MSCALAR,i,j,k) !DEV: should not matter - to be removed? - KW
#else
              call Eos_getAbarZbar(solnVec=solnData(:,i,j,k),Ye=Ye)           !DEV: for now
#endif
#endif
              ekinElecFrac = Ye * hy_eMassInUAmu
              ekinIonFrac = 1 - ekinElecFrac
              ! update the ion energy
              aux1 = - dtdx(j) * (e1flx2 - e1flx1)                             &
                   + dt*0.5e00 * (1 - eKinElecFrac)                            &
                   *( rho_o*vely_o*grav1d_o                                  &
                     +solnData(DENS_VAR,i,j,k)*solnData(VELY_VAR,i,j,k)*grav1d &
                    )

#ifdef E1_VAR
              etot = (rho_o * solnData(  E1_VAR,i,j,k) + aux1)*inv_new_dens
#endif
              
#ifdef EION_VAR
              ! get the ion internal energy
              einternal = solnData(EION_VAR,i,j,k)

              ! update ion internal energy
              aux1 = -dtdx(j) * ( (eint1flx2 - eint1flx1)               &
                     -0.5e0*(vely_o + solnData(VELY_VAR,i,j,k))        &
                     *aold_t * (p1av2 - p1av1) )
              
!              einternal = (rho_o * einternal + aux1)*inv_new_dens
              einternal = (rho_o * einternal + rescaledEIonAdvectedPlusPdVV)*inv_new_dens
!!$              einternal = max(einternal, hy_smallp*inv_new_dens)

! compute the new kinetic energy       
              ekin = 0.5e0 * ( solnData(VELX_VAR,i,j,k)**2             &
                              +solnData(VELY_VAR,i,j,k)**2             &
                              +solnData(VELZ_VAR,i,j,k)**2)
              ekin = ekin * (1 - ekinElecFrac)

! test whether we should use the internal energy from the evolution
              if (.TRUE. .OR. einternal .LT. hy_eint1Switch*ekin) then
                 etot = einternal + ekin
              else
                 einternal = etot - ekin
              endif
              
              if (inShock) then
                 eintIncreaseV = solnData(DENS_VAR,i,j,k)*einternal - rho_o*solnData(EION_VAR,i,j,k)
                 ! update internal energy advected
                 aux1 = -dtdx(j) * (eia1flx2 - eia1flx1)

                 einternalAdvectedV = aux1

                 eionIncreaseAboveAdvectedV =  eintIncreaseV - einternalAdvectedV
#ifdef DEBUG_XHYDRO
997              format(1x,I2,' inShock IonsY',3(1x,1PG22.15))
                 print 997,i,eintIncreaseV,einternalAdvectedV,eionIncreaseAboveAdvectedV
#endif
              end if

#ifdef E1_VAR
              solnData(E1_VAR,i,j,k) = etot
#endif
              solnData(EION_VAR,i,j,k) = einternal
#else
              
#ifdef E1_VAR
              solnData(E1_VAR,i,j,k) = etot
#endif
#endif
!********** E1_VAR END **********

!********** E2_VAR BEGIN **********
#ifdef FLASH_MULTISPECIES
              call Eos_getAbarZbar(solnVec=solnData(:,i,j,k),Ye=Ye)           !DEV: for now
#else
#if defined(SUMY_MSCALAR) && defined(YE_MSCALAR)
              Ye = solnData(YE_MSCALAR,i,j,k) !DEV: should not matter - to be removed? - KW
#else
              call Eos_getAbarZbar(solnVec=solnData(:,i,j,k),Ye=Ye)           !DEV: for now
#endif
#endif
              ekinElecFrac = Ye * hy_eMassInUAmu
              ! update the electron energy
              aux1 = - dtdx(j) * (e2flx2 - e2flx1)                             &
                   + dt*0.5e00 * ekinElecFrac                                &
                   *( rho_o*vely_o*grav1d_o                                  &
                     +solnData(DENS_VAR,i,j,k)*solnData(VELY_VAR,i,j,k)*grav1d &
                    )

#ifdef E2_VAR
              etot = (rho_o * solnData(  E2_VAR,i,j,k) + aux1)*inv_new_dens
#endif

#ifdef EELE_VAR
              ! get the electron internal energy
              einternal = solnData(EELE_VAR,i,j,k)

              ! update electron internal energy
              aux1 = -dtdx(j) * ( (eint2flx2 - eint2flx1)               &
                     -0.5e0*(vely_o + solnData(VELY_VAR,i,j,k))        &
                     *aold_t * (p2av2 - p2av1) )
              
!              einternal = (rho_o * einternal + aux1)*inv_new_dens
              einternal = (rho_o * einternal + rescaledEEleAdvectedPlusPdVV)*inv_new_dens
!!$              einternal = max(einternal, hy_smallp*inv_new_dens)

! compute the new kinetic energy       
              ekin = 0.5e0 * ( solnData(VELX_VAR,i,j,k)**2             &
                              +solnData(VELY_VAR,i,j,k)**2             &
                              +solnData(VELZ_VAR,i,j,k)**2)
              ekin = ekin * ekinElecFrac

! test whether we should use the internal energy from the evolution
              if (.TRUE. .OR. einternal .LT. hy_eint2Switch*ekin) then
                 etot = einternal + ekin
              else
                 einternal = etot - ekin
              endif
              
              if (inShock .AND. .FALSE.) then
                 eionTargetIncreaseV = max(eionIncreaseAboveAdvectedV, ekinIonFrac*eintIncreaseAboveAdvectedV)
                 eionTargetIncreaseV = min(eionTargetIncreaseV, einternal - hy_smallp)
!!$                 eeleTargetIncreaseV = max(eionIncreaseAboveAdvectedV, ekinEleFrac*eintIncreaseAboveAdvectedV)
                 eionAdjustment = max(0.0,eionTargetIncreaseV-eionIncreaseAboveAdvectedV) * inv_new_dens
                 if (eionAdjustment .LE. 0.0) then
#ifdef DEBUG_XHYDRO
                    print*,i,' InShock ignoring negative eionAdjustment:',eionAdjustment
#endif
                    eionAdjustment = 0.0
                    inShock = .FALSE.
                 else
                    solnData(EION_VAR,i,j,k) = solnData(EION_VAR,i,j,k) + eionAdjustment
#ifdef E1_VAR
                    solnData(E1_VAR,i,j,k) = solnData(E1_VAR,i,j,k) + eionAdjustment
#endif
                    eeleAdjustment = - min(eionAdjustment, einternal)
#ifdef DEBUG_XHYDRO
                    if (eionAdjustment + eeleAdjustment > 0.0) then
                       print*,i,' InShock WARNING some of adjustments:',eionAdjustment+eeleAdjustment
                    end if
#endif
                    einternal = einternal + eeleAdjustment
                    if (einternal < hy_smallp*inv_new_dens) then
                       print*,i,' InShock WARNING electron energy now very small:',einternal
                    end if
                    etot = etot + eeleAdjustment
                 end if
              end if

#ifdef E2_VAR
              solnData(E2_VAR,i,j,k) = etot
#endif
              solnData(EELE_VAR,i,j,k) = einternal
#else
              
#ifdef E2_VAR
              solnData(E2_VAR,i,j,k) = etot
#endif
#endif
!********** E2_VAR END **********

!********** E3_VAR BEGIN **********
              ! update the radiation energy
              aux1 = - dtdx(j) * (e3flx2 - e3flx1)                             &
                   + dt*0.5e00  * 0.0                                              &
                   *( rho_o*vely_o*grav1d_o                                  &
                     +solnData(DENS_VAR,i,j,k)*solnData(VELY_VAR,i,j,k)*grav1d &
                    )

#ifdef E3_VAR
              etot = (rho_o * solnData(  E3_VAR,i,j,k) + aux1)*inv_new_dens
#endif
              
#ifdef ERAD_VAR
              ! get the radiation internal energy
              einternal = solnData(ERAD_VAR,i,j,k)

              ! update radiation internal energy
              aux1 = -dtdx(j) * ( (eint3flx2 - eint3flx1)               &
                     -0.5e0*(vely_o + solnData(VELY_VAR,i,j,k))        &
                     *aold_t * (p3av2 - p3av1) )
              
!              einternal = (rho_o * einternal + aux1)*inv_new_dens
              einternal = (rho_o * einternal + rescaledERadAdvectedPlusPdVV)*inv_new_dens
!!$              einternal = max(einternal, hy_smallp*inv_new_dens)

! compute the new kinetic energy        - ZERO for now (and perhaps forever)
              ekin = 0.0

! test whether we should use the internal energy from the evolution
!!$              if (einternal .LT. hy_eint3Switch*ekin) then
                 etot = einternal + ekin
!!$              else
!!$                 einternal = etot - ekin
!!$              endif
              
#ifdef E3_VAR
              solnData(E3_VAR,i,j,k) = etot
#endif
              solnData(ERAD_VAR,i,j,k) = einternal
#else
              
#ifdef E3_VAR
              solnData(E3_VAR,i,j,k) = etot
#endif
#endif
!********** E3_VAR END **********
              
              call internal_shiftEints(solnData(EION_VAR,i,j,k), &
                                    solnData(EELE_VAR,i,j,k), &
                                    solnData(ERAD_VAR,i,j,k), &
                                    0.0,0.0,0.0, &
                                    1.,1.,1.,i,j,k)
           end do
        end do
     end do
     
  end if

!===============================================================================

! update in z direction
      
  if ((NDIM == 3) .and. (xyzswp == SWEEP_Z)) then

     select case (rangeSwitch)
     case (UPDATE_INTERIOR)
        kmin  = blkLimits(LOW,KAXIS) + 1 !NGUARD*K3D + 2
        kmax  = blkLimits(HIGH,KAXIS) - 1 ! NGUARD*K3D + NZB - 1
        kskip = 1
     case (UPDATE_BOUND)
        kmin  = blkLimits(LOW,KAXIS) ! NGUARD*K3D + 1
        kmax  = blkLimits(HIGH,KAXIS) ! NGUARD*K3D + NZB
        kskip = blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS) ! NZB - 1
     case default
        kmin  = blkLimits(LOW,KAXIS) ! NGUARD*K3D + 1
        kmax  = blkLimits(HIGH,KAXIS) ! NGUARD*K3D + NZB
        kskip = 1
     end select

     do k = kmin, kmax, kskip
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS) ! NGUARD*K2D+1, NGUARD*K2D+NYB
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS) ! NGUARD+1, NGUARD+NXB

              dtdx(k)  = tempDtDx(i,j,k)

              rhoflx1  = tempFlx(RHO_FLUX,i,j,k)
              rhoflx2  = tempFlx(RHO_FLUX,i,j,k+1)

              uflx1    = tempFlx(U_FLUX,i,j,k)
              uflx2    = tempFlx(U_FLUX,i,j,k+1)

              pav1     = tempFlx(P_FLUX,i,j,k)
              pav2     = tempFlx(P_FLUX,i,j,k+1)

              p1av1     = tempFlx(PION_FLUX,i,j,k)
              p1av2     = tempFlx(PION_FLUX,i,j,k+1)
              p2av1     = tempFlx(PELE_FLUX,i,j,k)
              p2av2     = tempFlx(PELE_FLUX,i,j,k+1)
              p3av1     = tempFlx(PRAD_FLUX,i,j,k)
              p3av2     = tempFlx(PRAD_FLUX,i,j,k+1)

              utflx1   = tempFlx(UT_FLUX,i,j,k)
              utflx2   = tempFlx(UT_FLUX,i,j,k+1)

              uttflx1  = tempFlx(UTT_FLUX,i,j,k)
              uttflx2  = tempFlx(UTT_FLUX,i,j,k+1)

              eflx1    = tempFlx(E_FLUX,i,j,k)
              eflx2    = tempFlx(E_FLUX,i,j,k+1)

              e1flx1    = tempFlx(E1_FLUX,i,j,k)
              e1flx2    = tempFlx(E1_FLUX,i,j,k+1)
              e2flx1    = tempFlx(E2_FLUX,i,j,k)
              e2flx2    = tempFlx(E2_FLUX,i,j,k+1)
              e3flx1    = tempFlx(E3_FLUX,i,j,k)
              e3flx2    = tempFlx(E3_FLUX,i,j,k+1)

              eintflx1 = tempFlx(EINT_FLUX,i,j,k)
              eintflx2 = tempFlx(EINT_FLUX,i,j,k+1)

              eint1flx1 = tempFlx(EION_FLUX,i,j,k)
              eint1flx2 = tempFlx(EION_FLUX,i,j,k+1)
              eint2flx1 = tempFlx(EELE_FLUX,i,j,k)
              eint2flx2 = tempFlx(EELE_FLUX,i,j,k+1)
              eint3flx1 = tempFlx(ERAD_FLUX,i,j,k)
              eint3flx2 = tempFlx(ERAD_FLUX,i,j,k+1)

              eiaflx1 = tempFlx(EIA_FLUX,i,j,k)
              eiaflx2 = tempFlx(EIA_FLUX,i,j,k+1)

              eia1flx1 = tempFlx(EI1A_FLUX,i,j,k)
              eia1flx2 = tempFlx(EI1A_FLUX,i,j,k+1)
              eia2flx1 = tempFlx(EI2A_FLUX,i,j,k)
              eia2flx2 = tempFlx(EI2A_FLUX,i,j,k+1)
              eia3flx1 = tempFlx(EI3A_FLUX,i,j,k)
              eia3flx2 = tempFlx(EI3A_FLUX,i,j,k+1)

              oneflx1 = tempFlx(ONE_FLUX,i,j,k)
              oneflx2 = tempFlx(ONE_FLUX,i,j,k+1)
#ifdef VOLD_FLUX
              voFlx(1) =  tempFlx(VOLD_FLUX,i,j,k)
              voFlx(2) = -tempFlx(VOLD_FLUX,i,j,k+1)
#endif

              do kk = 1,hy_numXn
                 xnflx1(kk) = tempFlx(SPECIES_FLUX_BEGIN+kk-1,i,j,k)
                 xnflx2(kk) = tempFlx(SPECIES_FLUX_BEGIN+kk-1,i,j,k+1)
              end do

              aold_t   = tempArea(i,j,k)
              grav1d_o = tempGrav1d_o(i,j,k)
              grav1d   = tempGrav1d(i,j,k)
              fict1d   = tempFict(i,j,k)

              rho_o    = solnData(DENS_VAR,i,j,k)
              prevMassV(0) = rho_o
              velz_o   = solnData(VELZ_VAR,i,j,k)

              if ( hy_useCmaAdvection .and. NSPECIES > 1  ) then

! update the partial mass densities and passive scalars * density

                 prevMassV(1) = 0.0; prevMassV(2) = 0.0
                 do n = 1, hy_numXn
                    solnData(SPECIES_BEGIN-1+n,i,j,k) =               &
                              rho_o*solnData(SPECIES_BEGIN-1+n,i,j,k) &
                            -dtdx(k)*(xnflx2(n) - xnflx1(n))
                    prevMassV(1) = prevMassV(1) + dtdx(k) * xnflx1(n)
                    prevMassV(2) = prevMassV(2) - dtdx(k) * xnflx2(n)
                 end do

! update the total mass density

                 solnData(DENS_VAR,i,j,k) = 0.e0

                 do n = 1, NSPECIES
                    solnData(DENS_VAR,i,j,k) =                        &
                    solnData(DENS_VAR,i,j,k) +                        &
                    solnData(SPECIES_BEGIN-1+n,i,j,k)
                 end do

! recover partial densities and passive scalars * density

                 inv_new_dens = 1.e0/solnData(DENS_VAR,i,j,k)

                 do n = 1, hy_numXn
                    solnData(SPECIES_BEGIN-1+n,i,j,k) =               &
                    solnData(SPECIES_BEGIN-1+n,i,j,k) * inv_new_dens
                 end do

              else

! update the density               
                 solnData(DENS_VAR,i,j,k) = solnData(DENS_VAR,i,j,k) - &
                                       dtdx(k) * (rhoflx2-rhoflx1)
                 prevMassV(1) =  dtdx(k) * rhoflx1
                 prevMassV(2) = -dtdx(k) * rhoflx2

! update the mass fractions and passive scalars
                 do n = 1, hy_numXn
                    solnData(SPECIES_BEGIN-1+n,i,j,k) =               &
                           ( rho_o*solnData(SPECIES_BEGIN-1+n,i,j,k)  &
                            -dtdx(k)*(xnflx2(n) - xnflx1(n))         &
                           )/solnData(DENS_VAR,i,j,k)
                 end do

              end if

! limit the density 
              solnData(DENS_VAR,i,j,k) =                              &
                   max(hy_smlrho, solnData(DENS_VAR,i,j,k))

              inv_new_dens = 1.e0/solnData(DENS_VAR,i,j,k)
!#ifdef CIP
!
!              if ( itrcr > 0 ) then
!
!                 ! CIP update and limit tracer
!
!                 do n = itrcr,itrcr
!                    solnData(n,i,j,k) = atan(solnData(n,i,j,k))/(0.9999d0*pi) + 0.5e0
!                    solnData(n,i,j,k) = max(0.e0, min(1.e0, solnData(n,i,j,k) ))
!                 end do
!
!              end if
!#endif

! update the velocities
              aux1 = -dtdx(k)*( (uflx2 - uflx1)                      &
                               +aold_t*(pav2 - pav1))                &
               +0.5e0*dt                                             &
               *( (solnData(DENS_VAR,i,j,k) + rho_o)*fict1d           &
                 +(rho_o*grav1d_o + solnData(DENS_VAR,i,j,k)*grav1d)  &
                )

              solnData(VELZ_VAR,i,j,k) =                              &
                         ( rho_o * solnData(VELZ_VAR,i,j,k)           &
                          +aux1                                      &
                         )*inv_new_dens
              
              solnData(VELX_VAR,i,j,k) =                              &
                         ( rho_o * solnData(VELX_VAR,i,j,k)           &
                          -dtdx(k) * (utflx2 - utflx1)               &
                         )*inv_new_dens
              
              solnData(VELY_VAR,i,j,k) =                              &
                         ( rho_o * solnData(VELY_VAR,i,j,k)           &
                          -dtdx(k) * (uttflx2 - uttflx1)             &
                         )*inv_new_dens

! update the total energy
              aux1 = - dtdx(k) * (eflx2 - eflx1)                             &
                   + dt*0.5e00                                               &
                   *( rho_o*velz_o*grav1d_o                                  &
                     +solnData(DENS_VAR,i,j,k)*solnData(VELZ_VAR,i,j,k)*grav1d &
                    )

              etot = (rho_o * solnData(ENER_VAR,i,j,k) + aux1)*inv_new_dens
               
#ifdef EINT_VAR
! get the internal energy
              einternal = solnData(EINT_VAR,i,j,k)
              prevEcompV(0,0) = rho_o * solnData(EINT_VAR,i,j,k)

! update internal energy
              aux1 = -dtdx(k) * ( (eintflx2 - eintflx1)               &
                     -0.5e0*(velz_o + solnData(VELZ_VAR,i,j,k))        &
                     *aold_t * (pav2 - pav1) )
               
              einternal = (rho_o * einternal + aux1)*inv_new_dens
              einternal = max(einternal, hy_smallp*inv_new_dens)
               
! compute the new kinetic energy
              ekin = 0.5e0 * ( solnData(VELX_VAR,i,j,k)**2             &
                              +solnData(VELY_VAR,i,j,k)**2             &
                              +solnData(VELZ_VAR,i,j,k)**2)
               
! test whether we should use the internal energy from the evolution
              if (einternal .LT. hy_eintSwitch*ekin) then
                 etot = einternal + ekin
              else
                 einternal = etot - ekin
              endif
              
#ifdef DBGS_VAR
              solndata(DBGS_VAR,i,j,k) = 0.5 * &
                   ( solndata(DBGS_VAR,i,j,k)  + tempFlx(SHOK_FLUX,i,j,k)+tempFlx(SHOK_FLUX,i,j,k+1) )
#endif
              inShock = ( tempFlx(SHOK_FLUX,i,j,k)+tempFlx(SHOK_FLUX,i,j,k+1) > 0 )
              if (hy_3Ttry_useShockDetect) then
                 forceRageLike = (.NOT. inShock)
              else
                 forceRageLike = .FALSE.
              end if

              eintIncreaseV = solnData(DENS_VAR,i,j,k)*einternal - rho_o*solnData(EINT_VAR,i,j,k)
              ! update internal energy advected
              aux1 = -dtdx(k) * (eiaflx2 - eiaflx1)
              einternalAdvectedV = aux1
              eintIncreaseAboveAdvectedV =  eintIncreaseV - einternalAdvectedV

              solnData(ENER_VAR,i,j,k) = etot
              solnData(EINT_VAR,i,j,k) = einternal
#else
              
              solnData(ENER_VAR,i,j,k) = etot
#endif

!!$              if(hy_3Ttry_B==0 .OR. hy_3Ttry_B==1 .OR. (hy_3Ttry_D .GE. 2.0 .AND. hy_3Ttry_E==1)) then
                 aux1 = -dtdx(k) * (eia3flx2 - eia3flx1)
                 eRadAdvectedV = aux1
                 aux1 = -dtdx(k) * (eia2flx2 - eia2flx1)
                 eEleAdvectedV = aux1
                 aux1 = -dtdx(k) * (eia1flx2 - eia1flx1)
                 eIonAdvectedV = aux1
!!$              else if(hy_3Ttry_B_rad==0 .OR. hy_3Ttry_B_rad==1) then
!!$                 aux1 = -dtdx(k) * (eia3flx2 - eia3flx1)
!!$                 eRadAdvectedV = aux1
!!$              end if

              if(hy_3Ttry_B==3 .OR. hy_3Ttry_B_rad==3) then
                 prevVolV(0) = 1.0
                 prevVolV(1) = dtdx(k) * voFlx(1)
                 prevVolV(2) = dtdx(k) * voFlx(2)
                 ! Note prevEcompV(0,0) was initialized above from solnData(EINT_VAR,i,j,k), before that changed.
                 prevEcompV(1,0) = rho_o * solnData(EION_VAR,i,j,k)
                 prevEcompV(2,0) = rho_o * solnData(EELE_VAR,i,j,k)
                 prevEcompV(3,0) = rho_o * solnData(ERAD_VAR,i,j,k)
                 do side = 1,2      !first left then right interface
                    select case (side)
                    case(1)
                       prevEcompV(0,side) = dtdx(k) * eiaflx1
                       prevEcompV(1,side) = dtdx(k) * eia1flx1
                       prevEcompV(2,side) = dtdx(k) * eia2flx1
                       prevEcompV(3,side) = dtdx(k) * eia3flx1
                    case(2)
                       prevEcompV(0,side) = - dtdx(k) * eiaflx2
                       prevEcompV(1,side) = - dtdx(k) * eia1flx2
                       prevEcompV(2,side) = - dtdx(k) * eia2flx2
                       prevEcompV(3,side) = - dtdx(k) * eia3flx2
                    end select
                    if (voFlx(side) .GE. 0.0) then ! Stuff is coming in on this side
                       ! ... keep it. 
                    else                           ! Stuff is going out on this side
                       prevMassV(0) = prevMassV(0) + prevMassV(side)
                       prevVolV(0)  = prevVolV(0)  + prevVolV(side)
                       prevEcompV(0:3,0)  = prevEcompV(0:3,0)  + prevEcompV(0:3,side)
                       prevMassV(side) = 0.0
                       prevVolV(side) = 0.0
                    end if
                 end do

                 newEinternalV(0:3) = 0.0
                 do piece = 0,2      !all fluid pieces contributing to new cell contents now
                    newVolV(piece) = prevMassV(piece) * inv_new_dens
                    if (newVolV(piece) .LE. 0.0) then
                       kappa(piece) = 1.0
                       prevEcompV(0:3,piece) = 0.0
                       newEcompV(0:3,piece) = 0.0
                    else
                       kappa(piece) = prevVolV(piece) / newVolV(piece) !compression factor
                       newEcompV(0:3,piece) = prevEcompV(0:3,piece) * kappa(piece)**(gamComp(0:3) - 1.0)
                       newEinternalV(0:3) = newEinternalV(0:3) + newEcompV(0:3,piece)
                    end if
                 end do
              end if

              if(hy_3Ttry_B==1) then 
                 aux1 = -dtdx(k) * (oneflx2 - oneflx1)
                 if(hy_3Ttry_G==1) then
                    eEleAdvectedPlusPdVV = eEleAdvectedV + 0.5*(p2av1+p2av2)*aux1
                    eIonAdvectedPlusPdVV = eIonAdvectedV + 0.5*(p1av1+p1av2)*aux1
                    eRadAdvectedPlusPdVV = eRadAdvectedV + 0.5*(p3av1+p3av2)*aux1
                 else
                    eEleAdvectedPlusPdVV = eEleAdvectedV + solnData(PELE_VAR,i,j,k)*aux1
                    eIonAdvectedPlusPdVV = eIonAdvectedV + solnData(PION_VAR,i,j,k)*aux1
                    eRadAdvectedPlusPdVV = eRadAdvectedV + solnData(PRAD_VAR,i,j,k)*aux1
                 end if
              else if(hy_3Ttry_B_rad==1) then 
                 aux1 = -dtdx(k) * (oneflx2 - oneflx1)
                 if(hy_3Ttry_G==1) then
                    eRadAdvectedPlusPdVV = eRadAdvectedV + 0.5*(p3av1+p3av2)*aux1
                 else
                    eRadAdvectedPlusPdVV = eRadAdvectedV + solnData(PRAD_VAR,i,j,k)*aux1
                 end if
              end if
              if(hy_3Ttry_B==2) then 
                 aux1 = -dtdx(k) * ( (eint1flx2 - eint1flx1)               &
                      -0.5e0*(velz_o + solnData(VELZ_VAR,i,j,k))        &
                      *aold_t * (p1av2 - p1av1) )
                 eIonAdvectedPlusPdVV =  aux1
                 aux1 = -dtdx(k) * ( (eint2flx2 - eint2flx1)               &
                      -0.5e0*(velz_o + solnData(VELZ_VAR,i,j,k))        &
                      *aold_t * (p2av2 - p2av1) )
                 eEleAdvectedPlusPdVV =  aux1
              end if
              if(hy_3Ttry_B_rad==2) then 
                 aux1 = -dtdx(k) * ( (eint3flx2 - eint3flx1)               &
                      -0.5e0*(velz_o + solnData(VELZ_VAR,i,j,k))        &
                      *aold_t * (p3av2 - p3av1) )
                 eRadAdvectedPlusPdVV =  aux1
              end if
!!#define yuckFactor (1.0)
              if(hy_3Ttry_B_rad==3) then
                 eRadAdvectedPlusPdVV =  newEinternalV(3) - rho_o*solnData(ERAD_VAR,i,j,k)
              end if
              if(hy_3Ttry_B_rad==0) then
                 rescaledERadAdvectedPlusPdVV = eRadAdvectedV
              else
                 rescaledERadAdvectedPlusPdVV = eRadAdvectedPlusPdVV
              endif
              if(hy_3Ttry_B==3) then 
                 eIonAdvectedPlusPdVV =  newEinternalV(1) - rho_o*solnData(EION_VAR,i,j,k)
                 eEleAdvectedPlusPdVV =  newEinternalV(2) - rho_o*solnData(EELE_VAR,i,j,k)
                 if ( eintIncreaseV - (eEleAdvectedPlusPdVV+rescaledERadAdvectedPlusPdVV)  &
                      < (eIonAdvectedV) + &
                      (eintIncreaseV - (eIonAdvectedV+eEleAdvectedV+eRadAdvectedV)) &
                      * (solnData(PION_VAR,i,j,k) / solnData(PRES_VAR,i,j,k)) &
                      ) then
                    forceRageLike = .TRUE.
                 else if (eIonAdvectedPlusPdVV*yuckFactor  &
                      + eEleAdvectedPlusPdVV+rescaledERadAdvectedPlusPdVV   &
                      > eintIncreaseV) then
                    forceRageLike = .TRUE.
                 else if ( eEleAdvectedPlusPdVV+rescaledERadAdvectedPlusPdVV   &
                      > eintIncreaseV) then
                    forceRageLike = .TRUE.
                 end if
              end if
              !!              print*,'eEleZ',eEleAdvectedV,eEleAdvectedPlusPdVV
              if((.not. forceRageLike) .AND.(hy_3Ttry_E==2 .OR. &
                   ((hy_3Ttry_D/=2.0).AND.(hy_3Ttry_D/=3.0)))) then !recalibrate based on energy ratio
                 if(hy_3Ttry_B==0) then
                    rescaledEEleAdvectedPlusPdVV = eEleAdvectedV
                    rescaledEIonAdvectedPlusPdVV = eIonAdvectedV
                 else
                    rescaledEEleAdvectedPlusPdVV = eEleAdvectedPlusPdVV
                    rescaledEIonAdvectedPlusPdVV = eIonAdvectedPlusPdVV
                 endif
                 zeroDummy = 0.0
                 if(hy_3Ttry_D/=0.0) then
                    if(hy_3Ttry_D==1.0) then
                       call Hydro_recalibrateEints(eintIncreaseV-rescaledEEleAdvectedPlusPdVV-rescaledERadAdvectedPlusPdVV,&
                            rescaledEIonAdvectedPlusPdVV,zeroDummy)
                    end if
                    if(hy_3Ttry_D==1.5) then
                       call Hydro_recalibrateEints(eintIncreaseV-rescaledEIonAdvectedPlusPdVV-rescaledERadAdvectedPlusPdVV,&
                            zeroDummy,rescaledEEleAdvectedPlusPdVV)
                    end if
                    if(hy_3Ttry_D==1.25) then
                       if (eintIncreaseV .GE. &
                            rescaledEEleAdvectedPlusPdVV+rescaledERadAdvectedPlusPdVV) then
                          call Hydro_recalibrateEints(eintIncreaseV-rescaledEEleAdvectedPlusPdVV-rescaledERadAdvectedPlusPdVV,&
                               rescaledEIonAdvectedPlusPdVV,zeroDummy)
                       else
                          call Hydro_recalibrateEints(eintIncreaseV,&
                               rescaledEIonAdvectedPlusPdVV,rescaledEEleAdvectedPlusPdVV,rescaledERadAdvectedPlusPdVV)
                       end if
                    end if
                    if(hy_3Ttry_D==1.75) then
                       call Hydro_recalibrateEints(eintIncreaseV-rescaledEEleAdvectedPlusPdVV,&
                            rescaledEIonAdvectedPlusPdVV,zeroDummy,rescaledERadAdvectedPlusPdVV)
                    end if
                    if(hy_3Ttry_D==1.875) then
                       call Hydro_recalibrateEints(eintIncreaseV-rescaledEIonAdvectedPlusPdVV,&
                            zeroDummy,rescaledEEleAdvectedPlusPdVV,rescaledERadAdvectedPlusPdVV)
                    end if
                    if(hy_3Ttry_D==2.0) then
                       if(hy_3Ttry_F==2) then
                          call Hydro_recalibrateEints(eintIncreaseV,&
                               rescaledEIonAdvectedPlusPdVV,rescaledEEleAdvectedPlusPdVV,rescaledERadAdvectedPlusPdVV)
                       end if
                    end if
                    if(hy_3Ttry_D==3.0) then
                       if(hy_3Ttry_F==2) then
                          call Hydro_recalibrateEints(eintIncreaseV-rescaledERadAdvectedPlusPdVV,&
                               rescaledEIonAdvectedPlusPdVV,rescaledEEleAdvectedPlusPdVV,zeroDummy)
                       end if
                    end if
                 end if
              else if(hy_3Ttry_E==1 .OR. forceRageLike) then  !recalibrate based on pressure ratio
                 zeroDummy = 0.0
                 PeP = solnData(PELE_VAR,i,j,k) / solnData(PRES_VAR,i,j,k)
                 PiP = solnData(PION_VAR,i,j,k) / solnData(PRES_VAR,i,j,k)
                 if (hy_3Ttry_D==3.0) then
                    PrP = 0.0
                 else
                    PrP = solnData(PRAD_VAR,i,j,k) / solnData(PRES_VAR,i,j,k)
                 end if
                 eintIncreaseAboveAdvectedV = eintIncreaseV - (eIonAdvectedV+eEleAdvectedV+eRadAdvectedV)
                 call Hydro_recalibrateEints(eintIncreaseAboveAdvectedV,&
                      PiP,PeP,PrP)
                 rescaledEEleAdvectedPlusPdVV = eEleAdvectedV + PeP
                 rescaledEIonAdvectedPlusPdVV = eIonAdvectedV + PiP
                 if (hy_3Ttry_D .NE. 3.0 .OR. hy_3Ttry_B_rad==0) then
                    rescaledERadAdvectedPlusPdVV = eRadAdvectedV + PrP
                 else
                    rescaledERadAdvectedPlusPdVV = eRadAdvectedPlusPdVV
                 end if
                 call internal_shiftEints(rescaledEIonAdvectedPlusPdVV, &
                                       rescaledEEleAdvectedPlusPdVV, &
                                       rescaledERadAdvectedPlusPdVV, &
                                       solnData(EION_VAR,i,j,k),solnData(EELE_VAR,i,j,k), &
                                       solnData(ERAD_VAR,i,j,k), &
                                       rho_o, solnData(DENS_VAR,i,j,k),inv_new_dens,i,j,k)
              end if
#ifdef XXHYDRO
              print 800,'eEleZ',eEleAdvectedV,eEleAdvectedPlusPdVV,rescaledEEleAdvectedPlusPdVV,&
                   eIonAdvectedV,eIonAdvectedPlusPdVV,rescaledEIonAdvectedPlusPdVV
#endif


!********** E1_VAR BEGIN **********
#ifdef FLASH_MULTISPECIES
              call Eos_getAbarZbar(solnVec=solnData(:,i,j,k),Ye=Ye)           !DEV: for now
#else
#if defined(SUMY_MSCALAR) && defined(YE_MSCALAR)
              Ye = solnData(YE_MSCALAR,i,j,k) !DEV: should not matter - to be removed? - KW
#else
              call Eos_getAbarZbar(solnVec=solnData(:,i,j,k),Ye=Ye)           !DEV: for now
#endif
#endif
              ekinElecFrac = Ye * hy_eMassInUAmu
              ekinIonFrac = 1 - ekinElecFrac
              ! update the ion energy
              aux1 = - dtdx(k) * (e1flx2 - e1flx1)                             &
                   + dt*0.5e00 * (1 - eKinElecFrac)                            &
                   *( rho_o*velz_o*grav1d_o                                  &
                     +solnData(DENS_VAR,i,j,k)*solnData(VELZ_VAR,i,j,k)*grav1d &
                    )

#ifdef E1_VAR
              etot = (rho_o * solnData(  E1_VAR,i,j,k) + aux1)*inv_new_dens
#endif
              
#ifdef EION_VAR
              ! get the ion internal energy
              einternal = solnData(EION_VAR,i,j,k)

              ! update ion internal energy
              aux1 = -dtdx(k) * ( (eint1flx2 - eint1flx1)               &
                     -0.5e0*(velz_o + solnData(VELZ_VAR,i,j,k))        &
                     *aold_t * (p1av2 - p1av1) )
              
!              einternal = (rho_o * einternal + aux1)*inv_new_dens
              einternal = (rho_o * einternal + rescaledEIonAdvectedPlusPdVV)*inv_new_dens
!!$              einternal = max(einternal, hy_smallp*inv_new_dens)

! compute the new kinetic energy       
              ekin = 0.5e0 * ( solnData(VELX_VAR,i,j,k)**2             &
                              +solnData(VELZ_VAR,i,j,k)**2             &
                              +solnData(VELZ_VAR,i,j,k)**2)
              ekin = ekin * (1 - ekinElecFrac)

! test whether we should use the internal energy from the evolution
              if (.TRUE. .OR. einternal .LT. hy_eint1Switch*ekin) then
                 etot = einternal + ekin
              else
                 einternal = etot - ekin
              endif
              
              if (inShock) then
                 eintIncreaseV = solnData(DENS_VAR,i,j,k)*einternal - rho_o*solnData(EION_VAR,i,j,k)
                 ! update internal energy advected
                 aux1 = -dtdx(k) * (eia1flx2 - eia1flx1)

                 einternalAdvectedV = aux1

                 eionIncreaseAboveAdvectedV =  eintIncreaseV - einternalAdvectedV
#ifdef DEBUG_XHYDRO
996              format(1x,I2,' inShock IonsZ',3(1x,1PG22.15))
                 print 996,i,eintIncreaseV,einternalAdvectedV,eionIncreaseAboveAdvectedV
#endif
              end if

#ifdef E1_VAR
              solnData(E1_VAR,i,j,k) = etot
#endif
              solnData(EION_VAR,i,j,k) = einternal
#else
              
#ifdef E1_VAR
              solnData(E1_VAR,i,j,k) = etot
#endif
#endif
!********** E1_VAR END **********

!********** E2_VAR BEGIN **********
#ifdef FLASH_MULTISPECIES
              call Eos_getAbarZbar(solnVec=solnData(:,i,j,k),Ye=Ye)           !DEV: for now
#else
#if defined(SUMY_MSCALAR) && defined(YE_MSCALAR)
              Ye = solnData(YE_MSCALAR,i,j,k) !DEV: should not matter - to be removed? - KW
#else
              call Eos_getAbarZbar(solnVec=solnData(:,i,j,k),Ye=Ye)           !DEV: for now
#endif
#endif
              ekinElecFrac = Ye * hy_eMassInUAmu
              ! update the electron energy
              aux1 = - dtdx(k) * (e2flx2 - e2flx1)                             &
                   + dt*0.5e00 * ekinElecFrac                                &
                   *( rho_o*velz_o*grav1d_o                                  &
                     +solnData(DENS_VAR,i,j,k)*solnData(VELZ_VAR,i,j,k)*grav1d &
                    )

#ifdef E2_VAR
              etot = (rho_o * solnData(  E2_VAR,i,j,k) + aux1)*inv_new_dens
#endif

#ifdef EELE_VAR
              ! get the electron internal energy
              einternal = solnData(EELE_VAR,i,j,k)

              ! update electron internal energy
              aux1 = -dtdx(k) * ( (eint2flx2 - eint2flx1)               &
                     -0.5e0*(velz_o + solnData(VELZ_VAR,i,j,k))        &
                     *aold_t * (p2av2 - p2av1) )
              
!              einternal = (rho_o * einternal + aux1)*inv_new_dens
              einternal = (rho_o * einternal + rescaledEEleAdvectedPlusPdVV)*inv_new_dens
!!$              einternal = max(einternal, hy_smallp*inv_new_dens)

! compute the new kinetic energy       
              ekin = 0.5e0 * ( solnData(VELX_VAR,i,j,k)**2             &
                              +solnData(VELY_VAR,i,j,k)**2             &
                              +solnData(VELZ_VAR,i,j,k)**2)
              ekin = ekin * ekinElecFrac

! test whether we should use the internal energy from the evolution
              if (.TRUE. .OR. einternal .LT. hy_eint2Switch*ekin) then
                 etot = einternal + ekin
              else
                 einternal = etot - ekin
              endif
              
              if (inShock .AND. .FALSE.) then
                 eionTargetIncreaseV = max(eionIncreaseAboveAdvectedV, ekinIonFrac*eintIncreaseAboveAdvectedV)
                 eionTargetIncreaseV = min(eionTargetIncreaseV, einternal - hy_smallp)
!!$                 eeleTargetIncreaseV = max(eionIncreaseAboveAdvectedV, ekinEleFrac*eintIncreaseAboveAdvectedV)
                 eionAdjustment = max(0.0,eionTargetIncreaseV-eionIncreaseAboveAdvectedV) * inv_new_dens
                 if (eionAdjustment .LE. 0.0) then
#ifdef DEBUG_XHYDRO
                    print*,i,' InShock ignoring negative eionAdjustment:',eionAdjustment
#endif
                    eionAdjustment = 0.0
                    inShock = .FALSE.
                 else
                    solnData(EION_VAR,i,j,k) = solnData(EION_VAR,i,j,k) + eionAdjustment
#ifdef E1_VAR
                    solnData(E1_VAR,i,j,k) = solnData(E1_VAR,i,j,k) + eionAdjustment
#endif
                    eeleAdjustment = - min(eionAdjustment, einternal)
#ifdef DEBUG_XHYDRO
                    if (eionAdjustment + eeleAdjustment > 0.0) then
                       print*,i,' InShock WARNING some of adjustments:',eionAdjustment+eeleAdjustment
                    end if
#endif
                    einternal = einternal + eeleAdjustment
                    if (einternal < hy_smallp*inv_new_dens) then
                       print*,i,' InShock WARNING electron energy now very small:',einternal
                    end if
                    etot = etot + eeleAdjustment
                 end if
              end if

#ifdef E2_VAR
              solnData(E2_VAR,i,j,k) = etot
#endif
              solnData(EELE_VAR,i,j,k) = einternal
#else
              
#ifdef E2_VAR
              solnData(E2_VAR,i,j,k) = etot
#endif
#endif
!********** E2_VAR END **********

!********** E3_VAR BEGIN **********
              ! update the radiation energy
              aux1 = - dtdx(k) * (e3flx2 - e3flx1)                             &
                   + dt*0.5e00  * 0.0                                              &
                   *( rho_o*velz_o*grav1d_o                                  &
                     +solnData(DENS_VAR,i,j,k)*solnData(VELZ_VAR,i,j,k)*grav1d &
                    )

#ifdef E3_VAR
              etot = (rho_o * solnData(  E3_VAR,i,j,k) + aux1)*inv_new_dens
#endif
              
#ifdef ERAD_VAR
              ! get the radiation internal energy
              einternal = solnData(ERAD_VAR,i,j,k)

              ! update radiation internal energy
              aux1 = -dtdx(k) * ( (eint3flx2 - eint3flx1)               &
                     -0.5e0*(velz_o + solnData(VELZ_VAR,i,j,k))        &
                     *aold_t * (p3av2 - p3av1) )
              
!              einternal = (rho_o * einternal + aux1)*inv_new_dens
              einternal = (rho_o * einternal + rescaledERadAdvectedPlusPdVV)*inv_new_dens
!!$              einternal = max(einternal, hy_smallp*inv_new_dens)

! compute the new kinetic energy        - ZERO for now (and perhaps forever)
              ekin = 0.0

! test whether we should use the internal energy from the evolution
!!$              if (einternal .LT. hy_eint3Switch*ekin) then
                 etot = einternal + ekin
!!$              else
!!$                 einternal = etot - ekin
!!$              endif
              
#ifdef E3_VAR
              solnData(E3_VAR,i,j,k) = etot
#endif
              solnData(ERAD_VAR,i,j,k) = einternal
#else
              
#ifdef E3_VAR
              solnData(E3_VAR,i,j,k) = etot
#endif
#endif
!********** E3_VAR END **********
              
              call internal_shiftEints(solnData(EION_VAR,i,j,k), &
                                    solnData(EELE_VAR,i,j,k), &
                                    solnData(ERAD_VAR,i,j,k), &
                                    0.0,0.0,0.0, &
                                    1.,1.,1.,i,j,k)
           end do
        end do
     end do
      
  end if

!=============================================================================

   return
 contains
   subroutine internal_shiftEints(deltaEion,deltaEele, &
                               deltaErad, &
                               eintIon_o,eintEle_o, &
                               eintRad_o, &
                               rho_o, newDens, inv_new_dens, i,j,k)

     use Hydro_data, ONLY : hy_smallEion,hy_smallEele,hy_smallErad

     real,intent(INOUT)          :: deltaEion,deltaEele
     real,intent(INOUT),OPTIONAL :: deltaErad
     real,intent(in),OPTIONAL    :: eintIon_o,eintEle_o,eintRad_o
     real,intent(in),OPTIONAL    :: rho_o, newDens, inv_new_dens
     integer,intent(in),OPTIONAL :: i,j,k

     real :: oldRho, newRho, newRhoInv
     real :: eIonNewAbove, eEleNewAbove, eRadNewAbove, eAllNewAbove
     real :: eIonUnder, eEleUnder, eRadUnder, eAllUnder
     if (present(rho_o)) then
        oldRho = rho_o
     else
        oldRho = 1.0
     end if
     if (present(newDens)) then
        newRho = newDens
     else
        newRho = oldRho
     end if
     if (present(inv_new_dens)) then
        newRhoInv = inv_new_dens
     else
        newRhoInv = 1.0 / newRho
     end if

     eIonNewAbove = (eintIon_o * oldRho + deltaEion)*newRhoInv - hy_smallEion
     eEleNewAbove = (eintEle_o * oldRho + deltaEele)*newRhoInv - hy_smallEele
     eRadNewAbove = (eintRad_o * oldRho + deltaErad)*newRhoInv - hy_smallErad
     eAllNewAbove = eIonNewAbove + eEleNewAbove + eRadNewAbove
     if (eAllNewAbove .LT. 0.0) then
        ! Cannot possibly shift energy around among components to remedy the
        ! combined energy deficit.
        return !RETURN, giving up!
     end if
     if (eIonNewAbove .GE. 0.0 .AND. eEleNewAbove .GE. 0.0) then 
        if (eRadNewAbove .GE. 0.0) then 
           ! Everything is peachy.
           return !RETURN, nothing to do!
        end if
     end if

#ifdef DEBUG_MAR2012
800  format(a,1P,3(G26.19),' @ (',i3,',',i3,',',i3,')')
     print 800,'_shift deltas BEFORE:',deltaEion,deltaEele,deltaErad,i,j,k
#endif

     eIonUnder = min(eintIon_o * oldRho + deltaEion - hy_smallEion*newRho, 0.0)
     eEleUnder = min(eintEle_o * oldRho + deltaEele - hy_smallEele*newRho, 0.0)
     eRadUnder = min(eintRad_o * oldRho + deltaErad - hy_smallErad*newRho, 0.0)

#ifdef DEBUG_MAR2012
     print 800,'_shift Unders BEFORE:',eIonUnder,eEleUnder,eRadUnder,i,j,k
#endif
     if (eIonUnder < 0.0) then
        deltaEion = deltaEion - eIonUnder
        if (eeleUnder == 0.0) then
#ifdef DEBUG_MAR2012
801        format('Shifting ',1P,G25.18,' from ',a,' to ',a,'...')
           print 801, -eIonUnder,'Electrons','Ions'
#endif
!!$           print '(a,2(G30.20))','deltaEele bef.:', deltaEele,eIonUnder
           deltaEele = deltaEele + eIonUnder !shift ion energy deficit to electrons
!!$           print '(a,2(G30.20))','deltaEele aft.:', deltaEele,spacing(deltaEele)
           eEleUnder = min(eintEle_o * oldRho + deltaEele - hy_smallEele*newRho, 0.0)
        else
#ifdef DEBUG_MAR2012
           print 801, -eIonUnder,'Radiation','Ions'
#endif
           deltaErad = deltaErad + eIonUnder !shift ion energy deficit to radiation
           eRadUnder = min(eintRad_o * oldRho + deltaErad - hy_smallErad*newRho, 0.0)
        end if
        eIonUnder = 0.0
     end if
     if (eEleUnder < 0.0) then  !last ditch effort:
#ifdef DEBUG_MAR2012
        print 801, -eEleUnder,'Radiation','Electrons'
#endif
        deltaEele = deltaEele - eEleUnder
        deltaErad = deltaErad + eEleUnder !shift electron energy deficit to radiation
        eRadUnder = min(eintRad_o * oldRho + deltaErad - hy_smallErad*newRho, 0.0)
        eEleUnder = 0.0
     end if


     if (eRadUnder < 0.0) then
        deltaErad = deltaErad - eRadUnder
        if (eeleUnder == 0.0) then
#ifdef DEBUG_MAR2012
           print 801, -eRadUnder,'Electrons','Radiation'
#endif
           deltaEele = deltaEele + eRadUnder !shift radiations energy deficit to electrons
!!$           eEleUnder = min(eintEle_o * oldRho + deltaEele - hy_smallEele*newRho, 0.0)
        else
#ifdef DEBUG_MAR2012
           print 801, -eRadUnder,'Ions','Radiation'
#endif
           deltaEion = deltaEion + eIonUnder !shift radiation energy deficit to ions
!!$           eIonUnder = min(eintIon_o * oldRho + deltaEion - hy_smallEion*newRho, 0.0)
        end if
     end if

#ifdef DEBUG_MAR2012
     print 800,'_shift deltas AFTER: ',deltaEion,deltaEele,deltaErad,i,j,k
#endif
   end subroutine internal_shiftEints

 end subroutine hy_ppm_updateSoln

