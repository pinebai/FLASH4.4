!!****if* source/physics/Hydro/HydroMain/split/PPM/threadBlockList/hy_ppm_block
!!
!! NAME
!!
!!  hy_ppm_block
!!
!! SYNOPSIS
!!
!!  call hy_ppm_block(integer(IN):: hy_meshMe,
!!                integer(IN):: sweepDir, 
!!                integer(IN):: blockID,
!!                real(IN)   :: dt, 
!!                real(IN)   :: dtOld, 
!!                integer(IN):: blkLimits(HIGH,MDIM),
!!                integer(IN):: blkLimitsGC(HIGH,MDIM),
!!                integer(IN):: bcs(2,MDIM),
!!                integer(IN):: numCells, 
!!                integer(IN):: numguard,
!!                real(IN)   :: primaryCoord(numCells),
!!                real(IN)   :: primaryLeftCoord(numCells),
!!                real(IN)   :: primaryRightCoord(numCells),
!!                real(IN)   :: primaryDx(numCells),
!!                real(IN)   :: secondCoord(numCells),
!!                real(IN)   :: thirdCoord(numCells),
!!                real(IN)   :: radialCoord(numCells),
!!                real(IN)   :: ugrid(numCells),
!!                real(OUT)  :: tempArea (blkLimitsGC(1,1):blkLimitsGC(2,1), blkLimitsGC(1,2):blkLimitsGC(2,2), blkLimitsGC(1,3):blkLimitsGC(2,3)),
!!                real(OUT)  :: tempGrav1d_o(blkLimitsGC(1,1):blkLimitsGC(2,1), blkLimitsGC(1,2):blkLimitsGC(2,2), blkLimitsGC(1,3):blkLimitsGC(2,3)),
!!                real(OUT)  :: tempGrav1d  (blkLimitsGC(1,1):blkLimitsGC(2,1), blkLimitsGC(1,2):blkLimitsGC(2,2), blkLimitsGC(1,3):blkLimitsGC(2,3)),
!!                real(OUT)  :: tempDtDx    (blkLimitsGC(1,1):blkLimitsGC(2,1), blkLimitsGC(1,2):blkLimitsGC(2,2), blkLimitsGC(1,3):blkLimitsGC(2,3)),
!!                real(OUT)  :: tempFict    (blkLimitsGC(1,1):blkLimitsGC(2,1), blkLimitsGC(1,2):blkLimitsGC(2,2), blkLimitsGC(1,3):blkLimitsGC(2,3)),
!!                real(OUT)  :: tempAreaLeft(blkLimitsGC(1,1):blkLimitsGC(2,1), blkLimitsGC(1,2):blkLimitsGC(2,2), blkLimitsGC(1,3):blkLimitsGC(2,3)),
!!                real(OUT)  :: tempFlx(NFLUXES, blkLimitsGC(1,1):blkLimitsGC(2,1), blkLimitsGC(1,2):blkLimitsGC(2,2), blkLimitsGC(1,3):blkLimitsGC(2,3)),
!!                real(IN)   :: shock       (blkLimitsGC(1,1):blkLimitsGC(2,1), blkLimitsGC(1,2):blkLimitsGC(2,2), blkLimitsGC(1,3):blkLimitsGC(2,3)),
!!                real, pointer :: solnData(:,:,:,:)
!!
!! DESCRIPTION
!!
!!  This routine takes a pointer to a block of grid data, picks out 1d
!!  slices, and applies the ppm algorithm, 1 row at a time, calculating
!!  fluxes.  The fluxes are stored and returned in the tempFlx/y/z
!!  arrays.  Other quantities used in the update step are calculated, as
!!  well, such as the geometry factors.
!!
!! ARGUMENTS
!!
!!   hy_meshMe -- my Processor Number
!!
!!   blockID -- My block number
!!
!!   sweepDir --  direction in which to do 1d hydro sweeps
!!
!!   dt --         timestep to advance through
!!
!!   dtOld --     previous timestep
!!
!!   blkLimits --  array holding upper and lower index limits of interior block cells (no GC)
!!
!!   blkLimitsGC --  array holding the upper and lower index limits of an entire block (including GC)
!!
!!   bcs --  boundary conditions for domain boundaries as obtained from Grid_getBlkBC
!!
!!   numCells --  tells this routine how long each 1d slice will be
!!
!!   numguard --  tells this routine how many guardcells surround the interior
!!
!!   primaryCoord --  positions of cell centers of the 1d slice
!!
!!   primaryLeftCoord --  positions of left interfaces
!!
!!   primaryRightCoord --  positions of right interfaces
!!
!!   primaryDx --  width of cells in the sweep direction
!!
!!   secondCoord --  for an x sweep: y coordinates; for a y sweep: x coord; for a z sweep: xcoord
!!
!!   thirdCoord --   for an x sweep: z coordinates; for a y sweep: z coord; for a z sweep: ycoord
!!
!!   radialCoord --  
!!
!!   ugrid --  this is for moving grid velocities
!!
!!   tempArea --  in case we're in non-cartesian geometry, compute cell face areas; used in the update step
!!
!!   tempGrav1d_o --  gravitational acceleration from the previous step, at cell edges
!!
!!   tempGrav1d --  gravitational acceleration from the current step 
!!
!!   tempDtDx --  dt/dx factor for updating quantities from fluxes
!!
!!   tempFict --  represents geometry related forces, e.g., centrifugal; used to update velocities 
!!
!!   tempFlx --   fluxes for all flux variables and all cells of the block; are used to update the solution, and
!!               given to Grid package for flux correction routines.
!!   tempAreaLeft --  cell face areas to the left
!!
!!
!!   shock --  0 if there is no shock or 1 if there is a shock for each cell; calculated by Hydro_detectShock 
!!
!!   solnData --  a pointer to the cell-centered data for the whole block
!!
!!***


! solnData depends on the ordering on unk
!!REORDER(4): solnData, tempFlx


#ifdef DEBUG_ALL
#define DEBUG_HYDRO
#endif


subroutine hy_ppm_block( hy_meshMe,blockID,sweepDir, dt, dtOld, &
                         blkLimits,blkLimitsGC,bcs,  &
                         numCells,numguard, &
                         primaryCoord ,     &
                         primaryLeftCoord , &
                         primaryRghtCoord , &
                         primaryDx        , &
                         secondCoord      , &
                         thirdCoord       , &
                         radialCoord      , &
                         ugrid            , &
                         tempArea, tempGrav1d_o, tempGrav1d,          &
                         tempDtDx, tempFict, tempAreaLeft,            &
                         tempFlx,                & 
                         shock, solnData )


  use Hydro_data, ONLY : hy_useGravity, hy_dirGeom,&
                         hy_dela, &
                         hy_dp , hy_du , hy_dut, hy_dutt, &
                         hy_drho, hy_dgame, hy_dgamc, hy_dgrav, &
                         hy_p6, hy_u6 , hy_ut6 , hy_utt6 , &
                         hy_rho6 , hy_game6 , hy_gamc6 , hy_grav6 , &
                         hy_pwl, hy_pwr, hy_dpw, hy_pw6l, hy_pw6r, hy_pwcubic,&
                         hy_gravl, hy_gravr,   &
                         hy_clft  , hy_plft  , hy_uttlft,   &
                         hy_ulft  , hy_vlft  , hy_utlft , &
                         hy_crght , hy_prght , hy_vrght ,         &
                         hy_urght , hy_utrght, hy_uttrgt, &
                         hy_gmelft, hy_gmergt, &
                         hy_gmclft, hy_gmcrgt, hy_pstor,hy_nriem, &
                         hy_deint,hy_eint6,hy_eiLft, hy_eiRght,     &
                         hy_dxn,hy_xn6,hy_xnlft, hy_xnrght,     &
                         hy_updateHydroFluxes, &
                         hy_useCellAreasForFluxes, &
                         hy_gravMass

  use Grid_interface, ONLY : Grid_getSingleCellVol, Grid_getBlkData
  use Driver_interface, ONLY : Driver_abortFlash
  use Timers_interface, ONLY : Timers_start, Timers_stop
  
  !! LOCAL - effective SCRATCH space
  !! shock and shock_multid are for shock detecting
  !! u -> utrt are for PPM.
  !! These are slices that we pass to PPM.  Optimized for vector
  !! machines years ago, PPM needs these slices in a very specific
  !! way.  This is costly.  

  implicit none
#include "constants.h"
#include "Flash.h"
#include "PPM.h"

  !! ------------
  !! ---- ARGUMENTS
  integer, intent(IN) :: hy_meshMe, blockID
  integer, intent(IN) :: sweepDir
  real,    intent(IN) :: dt, dtOld
  integer, intent(IN) :: numCells,numguard
  integer, intent(IN),dimension(2,MDIM) :: blkLimitsGC,blkLimits,bcs

  real,    pointer :: solnData(:,:,:,:) 

  !! ------------

  integer :: i1, i2, j1, j2, k1, k2
  integer :: i, j, k, ii,ilo,ihi,jlo,jhi,klo,khi
  integer :: iloGc,ihiGc,jloGc,jhiGc,kloGc,khiGc
  real    :: xbot, xtop, ybot, ytop, ylft, yrgt, zlft, zrgt
  integer :: numIntCells
  integer,parameter :: numXn=NSPECIES+NMASS_SCALARS
  integer :: sp,istat

  ! for gravity accumulation, implemented by LBR 12/19/2006
  integer, DIMENSION(MDIM) :: point
  real               :: cellVolume, tmpGravMass
  real, allocatable :: faceAreas(:,:,:), cellVolumes(:,:,:)

#ifdef FIXEDBLOCKSIZE
  real, intent(OUT), DIMENSION(GRID_ILO_GC:GRID_IHI_GC,       &
                               GRID_JLO_GC:GRID_JHI_GC,          &
                               GRID_KLO_GC:GRID_KHI_GC) ::       &
                               tempArea,       &
                               tempGrav1d_o,   &
                               tempGrav1d,     &
                               tempDtDx,       &
                               tempFict,       &
                               tempAreaLeft

  real, intent(IN), DIMENSION(GRID_ILO_GC:GRID_IHI_GC,       &
                              GRID_JLO_GC:GRID_JHI_GC,          &
                              GRID_KLO_GC:GRID_KHI_GC) :: &
                              shock

  real, intent(OUT), DIMENSION(NFLUXES,                   &
                               GRID_ILO_GC:GRID_IHI_GC,     &
                               GRID_JLO_GC:GRID_JHI_GC,     &
                               GRID_KLO_GC:GRID_KHI_GC) ::  &
                               tempFlx

  real,intent(IN), DIMENSION(MAXCELLS) :: primaryCoord ,  &
                                          primaryLeftCoord , &
                                          primaryRghtCoord , &
                                          primaryDx        , &
                                          secondCoord      , &
                                          thirdCoord       , &
                                          radialCoord     , &
                                          ugrid

  real, DIMENSION(MAXCELLS) :: dtdx, areaLeft, area, cvol, grav, ngrav, fict, shock_multid
  real, DIMENSION(MAXCELLS) :: rhoflx, uflx, pav, &
                               utflx, uttflx, eflx, eintflx,&
                               u, ut, utt, rho, p, e,  &
                               tmp, gamc, game,  &
                               uttp, utbt, utlt, utrt
  real, DIMENSION(MAXCELLS, NSPECIES+NMASS_SCALARS) :: xn,xnflx
  istat = 0
#else
  real, intent(OUT), DIMENSION(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),       &
                               blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),          &
                               blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) ::       &
                               tempArea,       &
                               tempGrav1d_o,   &
                               tempGrav1d,     &
                               tempDtDx,       &
                               tempFict,       &
                               tempAreaLeft
  real, intent(IN), DIMENSION(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
                              blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
                              blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) :: shock
  
  real, intent(OUT), DIMENSION(NFLUXES,                   &
                               blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),     &
                               blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),     &
                               blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) ::  &
                               tempFlx  
  real, intent(IN), DIMENSION(numCells) :: primaryCoord ,  &
                                           primaryLeftCoord , &
                                           primaryRghtCoord , &
                                           primaryDx        , &
                                           secondCoord      , &
                                           thirdCoord       , &
                                           radialCoord     , &
                                           ugrid

  real, DIMENSION(numCells) :: dtdx, areaLeft, area, cvol, grav, ngrav, fict, shock_multid
  real, DIMENSION(numCells) :: rhoflx, uflx, pav, &
                               utflx, uttflx, eflx, eintflx,&
                               u, ut, utt, rho, p, e,  &
                               tmp, gamc, game,  &
                               uttp, utbt, utlt, utrt
  real, DIMENSION(numCells, NSPECIES+NMASS_SCALARS) :: xn, xnflx

  allocate(hy_dela(numCells),stat = istat) 
  if (istat==0)allocate(hy_dp(numCells),stat = istat) 
  if (istat==0)allocate( hy_du(numCells),stat = istat)  
  if (istat==0)allocate(hy_dut(numCells),stat = istat) 
  if (istat==0)allocate( hy_dutt(numCells),stat = istat) 
  if (istat==0)allocate( hy_drho(numCells),stat = istat) 
  if (istat==0)allocate( hy_dgame(numCells),stat = istat) 
  if (istat==0)allocate( hy_dgamc(numCells),stat = istat) 
  if (istat==0)allocate( hy_dgrav(numCells),stat = istat) 
  if (istat==0)allocate( hy_p6(numCells),stat = istat) 
  if (istat==0)allocate( hy_u6(numCells),stat = istat)
  if (istat==0)allocate( hy_ut6(numCells),stat = istat) 
  if (istat==0)allocate( hy_utt6(numCells),stat = istat) 
  if (istat==0)allocate( hy_rho6(numCells),stat = istat) 
  if (istat==0)allocate( hy_game6(numCells),stat = istat) 
  if (istat==0)allocate( hy_gamc6(numCells),stat = istat) 
  if (istat==0)allocate( hy_grav6(numCells),stat = istat) 
  if (istat==0)allocate( hy_pwl(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_pwr(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_dpw(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_pw6l(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_pw6r(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_pwcubic(numCells),stat = istat) !! valid

  if (istat==0)allocate(hy_deint(numCells), stat = istat)
  if (istat==0)allocate(hy_eint6(numCells), stat = istat)
  if (istat==0)allocate(hy_eiLft(numCells),stat = istat) 
  if (istat==0)allocate( hy_eiRght(numCells),stat = istat) 

  if (istat==0)allocate(hy_dxn(numCells,numXN), stat = istat)
  if (istat==0)allocate(hy_xn6(numCells,numXN), stat = istat)
  if (istat==0)allocate(hy_xnlft(numCells,numXn),stat = istat) 
  if (istat==0)allocate( hy_xnrght(numCells,numXn),stat = istat) 

  if (istat==0)allocate(hy_gravl(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_gravr(numCells),stat = istat) !! valid
  if (istat==0)allocate(hy_clft(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_plft(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_uttlft(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_ulft(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_vlft(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_utlft(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_crght(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_prght(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_vrght(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_urght(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_utrght(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_uttrgt(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_gmelft(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_gmergt(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_gmclft(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_gmcrgt(numCells),stat = istat) !! valid
  
#endif

  if (istat==0)allocate( hy_pstor(hy_nriem + 2),stat = istat)
  if (istat .NE. 0) then
     call Driver_abortFlash("Memory allocation error in subroutine hy_block!")
  end if


!    Compute fluxes:  copy fluid data from 3D storage arrays        
!    into 1D work arrays, then operate on each row in turn.         
!    Update the solution in the 3D storage arrays, then save        
!    the computed fluxes for adjustment during the flux             
!    conservation step.                                             

!!$ IMPORTANT -- MOVE TO DATABASE. THIS IS A TEMPORARY HACK 
  grav(:) = 0.

  call Timers_start("hy_block")
! if we are using a hybrid Riemann solvers (i.e. using HLLE inside shocks), 
! then use a multi-dimensional shock detection -- this is more 
! accurate than the 1-d one done by the original PPM algorithm.
! Flags that signal in which cells a shock was detected are passed to this
! subroutine in the shock_multid array, if necessary.

  tempFlx(:,:,:,:) = 0.0
  select case (sweepDir)
  case (SWEEP_X)
     j1 = 1; j2 = 1
     k1 = 1; k2 = 1
     !! Loop over the interior to create the 1d slices
     iloGc = blkLimitsGC(LOW,IAXIS)
     ihiGc = blkLimitsGC(HIGH,IAXIS)
     ilo = blkLimits(LOW,IAXIS)
     ihi = blkLimits(HIGH,IAXIS)
     if (hy_useCellAreasForFluxes) then
        allocate(faceAreas(ihi+1, &
              blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
              blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)))
        call Grid_getBlkData(blockID, CELL_FACEAREA, ILO_FACE, EXTERIOR, &
                             (/1,blkLimits(LOW,JAXIS),blkLimits(LOW,KAXIS)/), &
                             faceAreas, &
          (/ihi+1, blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1, &
                   blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1/) )
        allocate(cellVolumes(ihi, &
              blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
              blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)))
        call Grid_getBlkData(blockID, CELL_VOLUME, 0, EXTERIOR, &
                             (/1,blkLimits(LOW,JAXIS),blkLimits(LOW,KAXIS)/), &
                             cellVolumes, &
          (/ihi, blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1, &
                 blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1/) )
     end if
     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        point(KAXIS) = k
        do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           point(JAXIS) = j
           if (hy_useCellAreasForFluxes) then
              areaLeft(1:ihi+1) = faceAreas(:,j,k)
              cvol(1:ihi) = cellVolumes(:,j,k)
           end if
           shock_multid(iloGc:ihiGc) = shock(iloGc:ihiGc,j,k)
           u(iloGc:ihiGc)    = solnData( VELX_VAR,iloGc:ihiGc, j, k )
           ut(iloGc:ihiGc)   = solnData( VELY_VAR, iloGc:ihiGc, j, k )
           utt(iloGc:ihiGc)  = solnData( VELZ_VAR, iloGc:ihiGc, j, k )
           rho(iloGc:ihiGc)  = solnData( DENS_VAR, iloGc:ihiGc, j, k )
           p(iloGc:ihiGc)    = solnData( PRES_VAR, iloGc:ihiGc, j, k )
           e(iloGc:ihiGc)    = solnData( ENER_VAR, iloGc:ihiGc, j, k )
           tmp(iloGc:ihiGc)  = solnData( TEMP_VAR, iloGc:ihiGc, j, k )
           game(iloGc:ihiGc) = solnData( GAME_VAR, iloGc:ihiGc, j, k )
           gamc(iloGc:ihiGc) = solnData( GAMC_VAR, iloGc:ihiGc, j, k )

           if (numXN > 0) then
#ifdef INDEXREORDER
              !! Note: solnData is going to be flipped so this is OK even if it does not seem so
              xn(iloGc:ihiGc,:) = solnData(SPECIES_BEGIN:(SPECIES_BEGIN+numXN-1),iloGc:ihiGc,j,k)
#else
              xn(iloGc:ihiGc,:) = transpose&
                (solnData(SPECIES_BEGIN:(SPECIES_BEGIN+numXN-1),iloGc:ihiGc,j,k ))
#endif
           end if
           if (NDIM.ge.2) then
              j1 = j - 1 ; j2 = j + 1
           end if
           if (NDIM.eq.3) then
              k1 = k - 1 ; k2 = k + 1
           end if

           !!Initialise to zero (because we may not be overwriting the entire array)
           uttp = 0
           utbt = 0
           utrt = 0
           utlt = 0
           uttp(iloGc:ihiGc) = solnData(VELY_VAR, iloGc:ihiGc, j2, k)
           utbt(iloGc:ihiGc) = solnData(VELY_VAR, iloGc:ihiGc, j1, k)
           utrt(iloGc:ihiGc) = solnData(VELZ_VAR, iloGc:ihiGc, j, k2)
           utlt(iloGc:ihiGc) = solnData(VELZ_VAR, iloGc:ihiGc, j, k1)

           xbot = 0.
           xtop = 0.           
           ybot = secondCoord(j1)
           ytop = secondCoord(j2)
           ylft = 0.
           yrgt = 0.
           zlft = thirdCoord(k1)
           zrgt = thirdCoord(k2)
           numIntCells = ihi-ilo+1

!        call Timers_start("hydro_1d")
           call hydro_1d (blockID,numIntCells,numCells,numguard, bcs,&
                          sweepDir, hy_meshMe,dt, dtOld, &
                          j, k,                       &
                          hy_dirGeom(IAXIS), hy_useGravity,              &
                          xbot, xtop,                 &
                          ybot, ytop, ylft, yrgt,     &
                          zlft, zrgt, ugrid,          &
                          primaryCoord ,              &
                          primaryLeftCoord ,          &
                          primaryRghtCoord ,          &
                          primaryDx        ,          &
                          secondCoord      ,          &
                          thirdCoord       ,          &
                          radialCoord     ,           &
                          u, ut, utt, rho, p, e, tmp, game, gamc,   &
                          xn, utbt, uttp, utlt, utrt,               &
                          shock_multid,                             &
                          dtdx, areaLeft, area, cvol, grav, ngrav, fict, &
                          rhoflx, uflx, pav, utflx, uttflx,         &
                          eflx, eintflx, xnflx)
!        call Timers_stop("hydro_1d")

           tempAreaLeft(ilo:ihi+1,j,k) = areaLeft(ilo:ihi+1)
         if (hy_updateHydroFluxes) then
           tmpGravMass = 0.0
           do i = ilo, ihi
              tempDtDx(i,j,k)     = dtdx(i)
              tempArea(i,j,k)     = area(i)
              tempGrav1d_o(i,j,k) = grav(i)
              tempGrav1d(i,j,k)   = ngrav(i)
              tempFict(i,j,k)     = fict(i)

              point(IAXIS) = i
              call Grid_getSingleCellVol(blockID,EXTERIOR,point,cellVolume)
              tmpGravMass = tmpGravMass + grav(i)*rho(i)*cellVolume
           enddo

           !$omp atomic
           hy_gravMass(IAXIS) = hy_gravMass(IAXIS) + tmpGravMass

           do i = ilo, ihi+1
              tempFlx(RHO_FLUX,i,j,k)  = rhoflx(i)
              tempFlx(U_FLUX,i,j,k)    = uflx(i)
              tempFlx(P_FLUX,i,j,k)    = pav(i)
              tempFlx(UT_FLUX,i,j,k)   = utflx(i)
              tempFlx(UTT_FLUX,i,j,k)  = uttflx(i)
              tempFlx(E_FLUX,i,j,k)    = eflx(i)
              tempFlx(EINT_FLUX,i,j,k) = eintflx(i)
           enddo
           
           do sp = 1, numXn
              do i = ilo, ihi+1
                 tempFlx(SPECIES_FLUX_BEGIN+sp-1,i,j,k) = xnflx(i,sp)
              end do
           enddo
         else
           do i = ilo, ihi
              tempDtDx(i,j,k)     = dtdx(i)
              tempArea(i,j,k)     = area(i)
              tempGrav1d_o(i,j,k) = grav(i)
              tempGrav1d(i,j,k)   = ngrav(i)
              tempFict(i,j,k)     = fict(i)
           end do
         end if
        
        end do !!j loop
     end do !!k loop

#if NDIM >= 2
  case (SWEEP_Y)
     k1 = 1; k2 = 1
     !! Loop over the interior to create the 1d slices
     jlo = blkLimits(LOW,JAXIS)
     jhi = blkLimits(HIGH,JAXIS)
     jloGc = blkLimitsGC(LOW,JAXIS)
     jhiGc = blkLimitsGC(HIGH,JAXIS)
     if (hy_useCellAreasForFluxes) then
        allocate(faceAreas(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                           jhi+1, &
                           blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)))
        call Grid_getBlkData(blockID, CELL_FACEAREA, JLO_FACE, EXTERIOR, &
                             (/blkLimits(LOW,IAXIS),1,blkLimits(LOW,KAXIS)/), &
                             faceAreas, &
                             (/blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1, &
                               jhi+1, &
                               blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1/) )
        allocate(cellVolumes(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                             jhi, &
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)))
        call Grid_getBlkData(blockID, CELL_VOLUME, 0, EXTERIOR, &
                             (/blkLimits(LOW,IAXIS),1,blkLimits(LOW,KAXIS)/), &
                             cellVolumes, &
                             (/blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1, &
                               jhi, &
                               blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1/) )
     end if
     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        point(KAXIS) = k
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           point(IAXIS) = i
           
           if (hy_useCellAreasForFluxes) then
              areaLeft(1:jhi+1) = faceAreas(i,:,k)
              cvol(1:jhi) = cellVolumes(i,:,k)
           end if
           shock_multid(jloGc:jhiGc) = shock(i,jloGc:jhiGc,k)
           u(jloGc:jhiGc)    = solnData( VELY_VAR,i, jloGc:jhiGc, k )
           ut(jloGc:jhiGc)   = solnData( VELX_VAR,i, jloGc:jhiGc, k )
           utt(jloGc:jhiGc)  = solnData( VELZ_VAR,i, jloGc:jhiGc, k )
           rho(jloGc:jhiGc)  = solnData( DENS_VAR,i, jloGc:jhiGc, k )
           p(jloGc:jhiGc)    = solnData( PRES_VAR,i, jloGc:jhiGc, k )
           e(jloGc:jhiGc)    = solnData( ENER_VAR,i, jloGc:jhiGc, k )
           tmp(jloGc:jhiGc)  = solnData( TEMP_VAR,i, jloGc:jhiGc, k )
           game(jloGc:jhiGc) = solnData( GAME_VAR,i, jloGc:jhiGc, k )
           gamc(jloGc:jhiGc) = solnData( GAMC_VAR,i, jloGc:jhiGc, k )
           if (numXN > 0) then
#ifdef INDEXREORDER
              xn(jloGc:jhiGc,:) = solnData(SPECIES_BEGIN:(SPECIES_BEGIN+numXN-1),i,jloGc:jhiGc,k)
#else
              xn(jloGc:jhiGc,:) = transpose&
                &(solnData(SPECIES_BEGIN:(SPECIES_BEGIN+numXN-1),i,jloGc:jhiGc,k ))
#endif
           end if
           i1 = i - 1 ; i2 = i + 1
           if (NDIM.eq.3) then
              k1 = k - 1 ; k2 = k + 1
           end if
           
           !!Initialise to zero (because we may not be overwriting the entire array)
           uttp = 0
           utbt = 0
           utrt = 0
           utlt = 0
           uttp(jloGc:jhiGc) = solnData(VELX_VAR, i2 ,jloGc:jhiGc, k)
           utbt(jloGc:jhiGc) = solnData(VELX_VAR, i1,jloGc:jhiGc, k)
           utrt(jloGc:jhiGc) = solnData(VELZ_VAR,i, jloGc:jhiGc, k2)
           utlt(jloGc:jhiGc) = solnData(VELZ_VAR,i, jloGc:jhiGc, k1)

           xbot = secondCoord(i1)
           xtop = secondCoord(i2)
           zlft = thirdCoord(k1)
           zrgt = thirdCoord(k2)
           numIntCells = jhi-jlo+1
!        call Timers_start("hydro_1d")
           call hydro_1d (blockID,numIntCells,numCells,numguard,bcs, &
                          sweepDir, hy_meshMe ,dt,dtOld, &
                          i, k,                       &
                          hy_dirGeom(JAXIS), hy_useGravity,                           &
                          xbot, xtop,                              &
                          ybot, ytop, ylft, yrgt,                  &
                          zlft, zrgt, ugrid,                       &
                          primaryCoord ,     &
                          primaryLeftCoord , &
                          primaryRghtCoord , &
                          primaryDx        , &
                          secondCoord      , &
                          thirdCoord       , &
                          radialCoord     , &
                          u, ut, utt, rho, p, e, tmp, game, gamc,  &
                          xn, utbt, uttp, utlt, utrt,              &
                          shock_multid,                            &
                          dtdx, areaLeft, area, cvol, grav, ngrav, fict, &
                          rhoflx, uflx, pav, utflx, uttflx,        &
                          eflx, eintflx, xnflx)
!        call Timers_stop("hydro_1d")

           tempAreaLeft(i,jlo:jhi+1,k) = areaLeft(jlo:jhi+1)
         if (hy_updateHydroFluxes) then
           tmpGravMass = 0.0
           do j = jlo, jhi
              tempDtDx(i,j,k)     = dtdx(j)
              tempArea(i,j,k)     = area(j)
              tempGrav1d_o(i,j,k) = grav(j)   !! Gradient in y direction
              tempGrav1d(i,j,k)   = ngrav(j)  !! Gradient in y direction + (dt/dtOld)*(difference in time)
              tempFict(i,j,k)     = fict(j)

              point(JAXIS) = j
              call Grid_getSingleCellVol(blockID,EXTERIOR,point,cellVolume)
              tmpGravMass = tmpGravMass + grav(j)*rho(j)*cellVolume !  grav(j)*rho(j)*cellVolume 

           enddo

           !$omp atomic
           hy_gravMass(JAXIS) = hy_gravMass(JAXIS) + tmpGravMass

           do j = jlo, jhi+1
              tempFlx(RHO_FLUX,i,j,k)  = rhoflx(j)
              tempFlx(U_FLUX,i,j,k)    = uflx(j)
              tempFlx(P_FLUX,i,j,k)    = pav(j)
              tempFlx(UT_FLUX,i,j,k)   = utflx(j)
              tempFlx(UTT_FLUX,i,j,k)  = uttflx(j)
              tempFlx(E_FLUX,i,j,k)    = eflx(j)
              tempFlx(EINT_FLUX,i,j,k) = eintflx(j)
           enddo
           
           do sp = 1, numXn
              do j = jlo, jhi+1
                 tempFlx(SPECIES_FLUX_BEGIN+sp-1,i,j,k) = xnflx(j,sp)
              end do
           enddo
         else
           do j = jlo, jhi
              tempDtDx(i,j,k)     = dtdx(j)
              tempArea(i,j,k)     = area(j)
              tempGrav1d_o(i,j,k) = grav(j)
              tempGrav1d(i,j,k)   = ngrav(j)
              tempFict(i,j,k)     = fict(j)
           end do
         end if
        end do !!j loop
     end do !!k loop

#endif
#if NDIM == 3 
  case (SWEEP_Z) 


     !! Loop over the interior to create the 1d slices
     klo = blkLimits(LOW,KAXIS)
     khi = blkLimits(HIGH,KAXIS)
     kloGc = blkLimitsGC(LOW,KAXIS)
     khiGc = blkLimitsGC(HIGH,KAXIS)

     if (hy_useCellAreasForFluxes) then
        allocate(faceAreas(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                           blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
                           khi+1))
        call Grid_getBlkData(blockID, CELL_FACEAREA, KLO_FACE, EXTERIOR, &
                             (/blkLimits(LOW,IAXIS),blkLimits(LOW,JAXIS),1/), &
                             faceAreas, &
                             (/blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1, &
                               blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1, &
                               khi+1/) )
        allocate(cellVolumes(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
                             khi))
        call Grid_getBlkData(blockID, CELL_VOLUME, 0, EXTERIOR, &
                             (/blkLimits(LOW,IAXIS),blkLimits(LOW,JAXIS),1/), &
                             cellVolumes, &
                             (/blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1, &
                               blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1, &
                               khi/) )
     end if
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        point(JAXIS) = j
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           point(IAXIS) = i
           
           if (hy_useCellAreasForFluxes) then
              areaLeft(1:khi+1) = faceAreas(i,j,:)
              cvol(1:khi) = cellVolumes(i,j,:)
           end if
           shock_multid(kloGc:khiGc) = shock(i,j,kloGc:khiGc)
           
           u(kloGc:khiGc)    = solnData( VELZ_VAR, i, j, kloGc:khiGc )
           ut(kloGc:khiGc)   = solnData( VELX_VAR, i, j, kloGc:khiGc )
           utt(kloGc:khiGc)  = solnData( VELY_VAR, i, j, kloGc:khiGc )
           rho(kloGc:khiGc)  = solnData( DENS_VAR, i, j, kloGc:khiGc )
           p(kloGc:khiGc)    = solnData( PRES_VAR, i, j, kloGc:khiGc )
           e(kloGc:khiGc)    = solnData( ENER_VAR, i, j, kloGc:khiGc )
           tmp(kloGc:khiGc)  = solnData( TEMP_VAR, i, j, kloGc:khiGc )
           game(kloGc:khiGc) = solnData( GAME_VAR, i, j, kloGc:khiGc )
           gamc(kloGc:khiGc) = solnData( GAMC_VAR, i, j, kloGc:khiGc )

           if (numXN > 0) then
#ifdef INDEXREORDER
              !! solnData will be flipped so this is OK
              xn(kloGc:khiGc,:) = solnData(SPECIES_BEGIN:(SPECIES_BEGIN+numXN-1), i, j, kloGc:khiGc)
#else
              xn(kloGc:khiGc,:) = transpose&
                &(solnData(SPECIES_BEGIN:(SPECIES_BEGIN+numXN-1), i, j, kloGc:khiGc))
#endif
           end if

           i1 = i - 1 ; i2 = i + 1
           j1 = j - 1 ; j2 = j + 1
           
           !!Initialise to zero (because we may not be overwriting the entire array)
           uttp = 0
           utbt = 0
           utrt = 0
           utlt = 0
           uttp(kloGc:khiGc) = solnData(VELX_VAR, i2, j, kloGc:khiGc)
           utbt(kloGc:khiGc) = solnData(VELX_VAR, i1, j, kloGc:khiGc)
           utrt(kloGc:khiGc) = solnData(VELY_VAR, i, j2, kloGc:khiGc)
           utlt(kloGc:khiGc) = solnData(VELY_VAR, i, j1, kloGc:khiGc)
           
           xbot = secondCoord(i1)
           xtop = secondCoord(i2)
           ylft = thirdCoord(j1)
           yrgt = thirdCoord(j2)
           numIntCells = khi-klo+1
!        call Timers_start("hydro_1d")
           call hydro_1d (blockID,numIntCells,numCells,numguard,bcs, &
                          sweepDir, hy_meshMe ,dt,dtOld, &
                          i, j,                       &
                          hy_dirGeom(KAXIS), hy_useGravity,                           &
                          xbot, xtop,                              &
                          ybot, ytop, ylft, yrgt,                  &
                          zlft, zrgt, ugrid,                       &
                          primaryCoord ,     &
                          primaryLeftCoord , &
                          primaryRghtCoord , &
                          primaryDx        , &
                          secondCoord      , &
                          thirdCoord       , &
                          radialCoord     , &
                          u, ut, utt, rho, p, e, tmp, game, gamc,  &
                          xn, utbt, uttp, utlt, utrt,              &
                          shock_multid,                            &
                          dtdx, areaLeft, area, cvol, grav, ngrav, fict, &
                          rhoflx, uflx, pav, utflx, uttflx,        &
                          eflx, eintflx, xnflx)
!        call Timers_stop("hydro_1d")
           tempAreaLeft(i,j,klo:khi+1) = areaLeft(klo:khi+1)
         if (hy_updateHydroFluxes) then
           tmpGravMass = 0.0
           do k = klo, khi 
              tempDtDx(i,j,k)     = dtdx(k)
              tempArea(i,j,k)     = area(k)
              tempGrav1d_o(i,j,k) = grav(k)
              tempGrav1d(i,j,k)   = ngrav(k)
              tempFict(i,j,k)     = fict(k)

              point(KAXIS) = k
              call Grid_getSingleCellVol(blockID,EXTERIOR,point,cellVolume)
              tmpGravMass = tmpGravMass + grav(k)*rho(k)*cellVolume
           enddo

           !$omp atomic
           hy_gravMass(KAXIS) = hy_gravMass(KAXIS) + tmpGravMass

           do k = klo, khi+1 
              tempFlx(RHO_FLUX,i,j,k)  = rhoflx(k)
              tempFlx(U_FLUX,i,j,k)    = uflx(k)
              tempFlx(P_FLUX,i,j,k)    = pav(k)
              tempFlx(UT_FLUX,i,j,k)  = utflx(k)
              tempFlx(UTT_FLUX,i,j,k)  = uttflx(k)
              tempFlx(E_FLUX,i,j,k)    = eflx(k)
              tempFlx(EINT_FLUX,i,j,k) = eintflx(k)
           enddo
           
           do sp = 1, numXn
              do k = klo, khi+1 
                 tempFlx(SPECIES_FLUX_BEGIN+sp-1,i,j,k) = xnflx(k,sp)
              end do
           enddo
           
         else
           do k = klo, khi
              tempDtDx(i,j,k)     = dtdx(k)
              tempArea(i,j,k)     = area(k)
              tempGrav1d_o(i,j,k) = grav(k)
              tempGrav1d(i,j,k)   = ngrav(k)
              tempFict(i,j,k)     = fict(k)
           end do
         end if
           
        end do
     end do
     
#endif
     
  end select

  if (hy_useCellAreasForFluxes) then
     deallocate(cellVolumes)
     deallocate(faceAreas)
  end if
  deallocate( hy_pstor)

#ifndef FIXEDBLOCKSIZE
  deallocate(hy_dela) 
  deallocate(hy_dp) 
  deallocate( hy_du)  
  deallocate(hy_dut) 
  deallocate( hy_dutt) 
  deallocate( hy_drho) 
  deallocate( hy_dgame) 
  deallocate( hy_dgamc) 
  deallocate( hy_dgrav) 
  deallocate( hy_p6) 
  deallocate( hy_u6)
  deallocate( hy_ut6) 
  deallocate( hy_utt6) 
  deallocate( hy_rho6) 
  deallocate( hy_game6) 
  deallocate( hy_gamc6) 
  deallocate( hy_grav6) 
  deallocate( hy_pwl) 
  deallocate( hy_pwr) 
  deallocate( hy_dpw) 
  deallocate( hy_pw6l) 
  deallocate( hy_pw6r) 
  deallocate( hy_pwcubic) 

  deallocate(hy_deint)
  deallocate(hy_eint6)
  deallocate(hy_eiLft)
  deallocate( hy_eiRght)

  deallocate(hy_dxn)
  deallocate(hy_xn6)
  deallocate(hy_xnlft)
  deallocate( hy_xnrght)

  deallocate(hy_gravl) 
  deallocate( hy_gravr) 
  deallocate(hy_clft) 
  deallocate( hy_plft) 
  deallocate( hy_uttlft) 
  deallocate( hy_ulft) 
  deallocate( hy_vlft) 
  deallocate( hy_utlft) 
  deallocate( hy_crght) 
  deallocate( hy_prght) 
  deallocate( hy_vrght) 
  deallocate( hy_urght) 
  deallocate( hy_utrght) 
  deallocate( hy_uttrgt) 
  deallocate( hy_gmelft) 
  deallocate( hy_gmergt) 
  deallocate( hy_gmclft) 
  deallocate( hy_gmcrgt) 
  
#endif
  call Timers_stop("hy_block")
end subroutine hy_ppm_block

