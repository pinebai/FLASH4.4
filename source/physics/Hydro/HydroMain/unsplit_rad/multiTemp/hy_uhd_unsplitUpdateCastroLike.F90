!!****if* source/physics/Hydro/HydroMain/unsplit_rad/multiTemp/hy_uhd_unsplitUpdateCastroLike
!!
!! NAME
!!
!!  hy_uhd_unsplitUpdateCastroLike
!!
!! SYNOPSIS
!!
!!  call hy_uhd_unsplitUpdateCastroLike( integer(IN) :: blockID,
!!                                 integer(IN) :: rangeSwitch,
!!                                 integer(IN) :: blkLimits(2,MDIM),
!!                                 integer(IN) :: dataSize(3),
!!                                 real(IN)    :: dt,
!!                                 real(IN)    :: del(MDM),
!!                                 real(IN)    :: xflux(:,:,:,:), 
!!                                 real(IN)    :: yflux(:,:,:,:), 
!!                                 real(IN)    :: zflux(:,:,:,:),
!!                                 real(IN)    :: gravX(3,:,:,:),
!!                                 real(IN)    :: gravY(3,:,:,:),
!!                                 real(IN)    :: gravZ(3,:,:,:))
!!
!! ARGUMENTS
!!
!!   blockID      - current block ID
!!   dt           - timestep
!!   del          - deltas in {x,y,z} directions
!!   dataSize     - size of the current block
!!   blkLimits    - an array that holds the lower and upper indices of the section
!!                  of block without the guard cells
!!   xflux,yflux,zflux - cell face centered fluxes at each {=x,y,z} direction
!!   gravX,gravY,gravZ - gravity components in x,y,z directions
!!
!! DESCRIPTION
!!
!!   This routine updates the cell-centered conservative variables and intermediate
!!   internal energy to the next time step using an directionally unsplit scheme.
!!
!! SIDE EFFECTS
!!
!!   This routine updates the following cell-centered conservative variables:
!!     EION_VAR
!!     EELE_VAR
!!     ERAD_VAR
!!
!!***

!!REORDER(4): U, scrch_Ptr, [xyz]flux

Subroutine hy_uhd_unsplitUpdateCastroLike &
     (blockID,rangeSwitch,blkLimits,datasize,dt,del,xflux,yflux,zflux, scrch_Ptr)

  use Hydro_data,           ONLY : hy_geometry
  use Hydro_data,           ONLY : hy_3TMode

  use Hydro_interface,      ONLY : Hydro_recalibrateEints

  use Driver_interface,     ONLY : Driver_abortFlash

  use hy_uhd_MultiTempData, ONLY : hy_3Ttry_B, &
                                   hy_3Ttry_E, &
                                   hy_3Ttry_F, &
                                   hy_3Ttry_G, &
                                   hy_3Ttry_D, &
                                   hy_3Ttry_B_rad, &
                                   hy_mtPresRatLambda3Min
  use hy_uhd_MultiTempData, ONLY : hy_mtScaleAccel, hy_mtScaleWork, hy_mtScaleLorentz
  use hy_uhd_MultiTempData, ONLY: fB0=>hy_mtFactorB0,fB1=>hy_mtFactorB1,fB2=>hy_mtFactorB2
  use hy_uhd_MultiTempData, ONLY: fB0r=>hy_mtFactorBrad0,fB1r=>hy_mtFactorBrad1,fB2r=>hy_mtFactorBrad2

  use Grid_interface,       ONLY : Grid_getBlkPtr, &
                                   Grid_releaseBlkPtr, &
                                   Grid_getBlkData, &
                                   Grid_getCellCoords, &
                                   Grid_getBlkIndexLimits
  
  use Logfile_interface,    ONLY : Logfile_stampMessage
  
  use Opacity_interface,    ONLY : Opacity

  use hy_uhd_interface,     ONLY : hy_uhd_ragelike

#include "constants.h"
#include "Flash.h"
#include "Eos.h"
#include "Eos_components.h"
#include "UHD.h"

#ifdef FLASH_USM_MHD
 use Hydro_data,           ONLY : hy_hallVelocity, hy_useMagneticResistivity
#ifdef FLASH_UHD_3T
  use hy_uhd_interface,ONLY : hy_uhd_getCurrents
#endif
#endif
#ifdef FLASH_UHD_3TRADFLAH
  use Hydro_data,          ONLY : hy_lam3ScaleFactor
#endif

 implicit none


  !! ---- Arguments ---------------------------------
  integer,intent(IN) :: blockID, rangeSwitch
  integer,intent(IN) :: blkLimits  (LOW:HIGH,MDIM)
  integer,intent(IN) :: datasize(MDIM)      
  real, intent(IN) :: dt
  real, intent(IN) :: del(MDIM)
#ifdef FIXEDBLOCKSIZE
  real, intent(in) :: xflux(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
  real, intent(in) :: yflux(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
  real, intent(in) :: zflux(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
#else
  real, intent(in) :: xflux(NFLUXES,datasize(IAXIS),datasize(JAXIS),datasize(KAXIS))  
  real, intent(in) :: yflux(NFLUXES,datasize(IAXIS),datasize(JAXIS),datasize(KAXIS))  
  real, intent(in) :: zflux(NFLUXES,datasize(IAXIS),datasize(JAXIS),datasize(KAXIS))
#endif
#if EOSCOMP_MATTER > 0
  logical,parameter :: do2T = .TRUE.  !!DEV: should check that it is 2.
#else
  logical,parameter :: do2T = .FALSE. !full 3T then!
#endif
  !!---------------------------------------------------

  real :: aux1,PeP, PiP, PrP, velx_left, velx_right, &
          dutot, &
          duele_adv,&
          duion_adv,&
          durad_adv,&
          duele, &
          duion, & 
          durad, &
          uion_adv, &
          uele_adv, &
          urad_adv, &
          eion_old, uion_old, uion_new, uion_ragelike, &
          eele_old, uele_old, uele_new, uele_ragelike, &
          erad_old, urad_old, urad_new, urad_ragelike, &
          dens_old, dens_new, densNewInv,&
          utot_old, utot_new, dutot_adv, &
          work_ele, work_rad, divv, Qohm
  real :: dueleAdvPlus, duionAdvPlus, duradAdvPlus, dutotAdvPlus
  real :: uref, duDirect
  real :: scaleWork, scaleAccel, scaleLorentz

  integer :: blklmts(LOW:HIGH,MDIM)
  integer :: blklmtsgc(LOW:HIGH,MDIM)
  
  integer :: i,j,k,iskip,kx,ky,kz
  integer :: imin,imax,jmin,jmax,kmin,kmax
  real    :: dx, dy, dz
  real, pointer, dimension(:,:,:,:) :: U
  real, pointer, dimension(:,:,:,:) :: scrch_Ptr ! this holds old density (VAR1) and old internal energy (VAR2)
  real, dimension(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)) :: faceAreas, cellVolumes
  real, dimension(3, dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)) :: Jp, Jm
  real, allocatable :: xcent(:), ycent(:), zcent(:)
  character(len=MAX_STRING_LENGTH) errmsg

  real :: lf, rf, lfj, rfj
  integer :: isize, jsize, ksize

  real    :: pele, pele_adv
  real    :: pion, pion_adv
  real    :: prad, prad_adv
  
  real :: kappaRatio
  real :: kappaP, kappaR, kappaDummy
  real :: lam3Rad, lam3RadFloored
  real :: uDotGradPradPhys, duradAccel, dueleAccel, duionAccel
  real :: lam3uDotGradPradPhys
  real :: duradWorkLike, dueleWorkLike, duionWorkLike
  real :: uDotGradPmat, duradLamPlus ,scaledDuradLamPlus
  real :: duradWork, dueleWork, duionWork, duradWorkPlus
  real :: presFL,presFR, presGL,presGR, presHL,presHR
  real :: pmatFL,pmatFR, pmatGL,pmatGR, pmatHL,pmatHR
  real :: pradFL,pradFR, pradGL,pradGR, pradHL,pradHR
  real :: pradEffFL,pradEffFR, pradEffGL,pradEffGR, pradEffHL,pradEffHR
  real :: velxFL,velxFR, velyGL,velyGR, velzHL,velzHR
  real :: lam3L,lam3R
  real :: lam3Lf,lam3Rf         !floored
  logical :: inShock

  if (hy_3TMode .NE. HY3T_CASTROLIKE) then
     call Driver_abortFlash('Dazed and confused')
  end if

  scaleWork = hy_mtScaleWork
  scaleAccel  = hy_mtScaleAccel
  scaleLorentz  = hy_mtScaleLorentz
  pele_adv = -1.0
  pion_adv = -1.0
  prad_adv = -1.0
  pele     = -1.0
  prad     = -1.0
  work_ele = -1.0
  work_rad = -1.0


  iSize = blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1
  jSize = blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1
  kSize = blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1
  
  Jp = 0.0 
  Jm = 0.0
  
  faceAreas = 0.
  call Grid_getBlkData(blockID, CELL_FACEAREA, ILO_FACE, EXTERIOR, &
       (/blkLimits(LOW,IAXIS),blkLimits(LOW,JAXIS),blkLimits(LOW,KAXIS)/), &
       faceAreas(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+1,&
                 blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),  &
                 blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)), &
       (/isize+1, jsize, ksize/) )

  cellVolumes = 0.
  call Grid_getBlkData(blockID, CELL_VOLUME, 0, EXTERIOR, &
       (/blkLimits(LOW,IAXIS),blkLimits(LOW,JAXIS),blkLimits(LOW,KAXIS)/), &
       cellVolumes(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                   blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
                   blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)), &
       (/isize, jsize, ksize/) )
  

  !! Set ranges for update
  imin  = blkLimits(LOW, IAXIS)
  imax  = blkLimits(HIGH,IAXIS)
  jmin  = 1
  jmax  = 1
  kmin  = 1
  kmax  = 1

  dx = del(DIR_X)
  dy = 1.
  dz = 1.
  if (NDIM >= 2) then
     jmin  = blkLimits(LOW, JAXIS)
     jmax  = blkLimits(HIGH,JAXIS)
     dy = del(DIR_Y)
     if (NDIM == 3) then
        kmin  = blkLimits(LOW, KAXIS)
        kmax  = blkLimits(HIGH,KAXIS)
        dz = del(DIR_Z)
     endif
  endif

  !! Set regions to update depending on update mode
  iskip = 1
  if (rangeSwitch==UPDATE_INTERIOR) then
     imin  = imin+1
     imax  = imax-1
     !iskip = 1
     if (NDIM >= 2) then
        jmin  = jmin+1
        jmax  = jmax-1
        if (NDIM == 3) then
           kmin  = kmin+1
           kmax  = kmax-1
        endif
     endif
  endif

  call Grid_getBlkIndexLimits(blockId,blklmts,blklmtsgc)
  
  allocate(xcent(blklmts(LOW, IAXIS) : blklmts(HIGH, IAXIS)))
  allocate(ycent(blklmts(LOW, JAXIS) : blklmts(HIGH, JAXIS)))
  allocate(zcent(blklmts(LOW, KAXIS) : blklmts(HIGH, KAXIS)))
  
  call Grid_getCellCoords(IAXIS, blockId, CENTER, .FALSE., xcent, size(xcent))
  call Grid_getCellCoords(JAXIS, blockId, CENTER, .FALSE., ycent, size(ycent))
  call Grid_getCellCoords(KAXIS, blockId, CENTER, .FALSE., zcent, size(zcent))
  

  ! Get block pointers
  call Grid_getBlkPtr(blockID,U,CENTER)

#ifdef FLASH_USM_MHD
#ifdef FLASH_UHD_3T
  if (hy_hallVelocity) then
  ! Get currents
    call hy_uhd_getCurrents(blockID, rangeSwitch, blkLimits,datasize, del, Jp, Jm, 1, &
                            scrch_Ptr)
  endif 
#endif
#endif

  ! define dimension dependent switches
  kx=0
  ky=0
  kz=0

  if (NDIM > 1) then
     ky=1
     if (NDIM > 2) then
        kz=1
     endif
  endif

  iskip = 1
  if (NDIM == 1 .and. rangeSwitch .eq. UPDATE_BOUND) iskip = imax-imin

  do k=kmin-kx*kz,kmax+kx*kz
     do j=jmin-kx*ky,jmax+kx*ky
        if (NDIM >= 2) then
           iskip = 1
           if (rangeSwitch == UPDATE_BOUND .and. j > jmin .and. j < jmax) then
              iskip = imax-imin
              if (NDIM == 3) then
                 iskip = 1
                 if (k > kmin .and. k < kmax) then
                    iskip = imax-imin
                 endif
              endif
           endif
        endif

        do i=imin-kx,imax+kx,iskip
           ! Loop over cells in block...

#ifdef BDRY_VAR
           if(U(BDRY_VAR,i,j,k) > 0.0) cycle
#endif

           ! **************************************************************
           ! *   Load variables into temporary variables for easier use   *
           ! **************************************************************
           dens_old = scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,i,j,k)
           dens_new = U(DENS_VAR,i,j,k)
           densNewInv = 1. / dens_new

           utot_new = U(DENS_VAR,i,j,k)*U(EINT_VAR,i,j,k)
           utot_old = scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,i,j,k)*scrch_Ptr(HY_VAR2_SCRATCHCTR_VAR,i,j,k)

           eele_old = U(EELE_VAR,i,j,k)
           eion_old = U(EION_VAR,i,j,k)
           erad_old = U(ERAD_VAR,i,j,k)

           uele_old = eele_old * dens_old
           uion_old = eion_old * dens_old
           urad_old = erad_old * dens_old
!!$print*,i,j,k
!!$!print*,dens_old,dens_new,utot_new,utot_old,eele_old,eion_old,erad_old
!!$print*,scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,i,j,k),scrch_Ptr(HY_VAR2_SCRATCHCTR_VAR,i,j,k)
!!$if (i==5) stop
           ! ******************************************
           ! *   Compute advected internal energies   *
           ! ******************************************

           ! Compute lf/rf -> these are geometric factors which will
           ! be 1 in Cartesian coordinates, but are needed in
           ! cylindrical coordinates...
           if(hy_geometry == CARTESIAN) then
              lf = 1.0
              rf = 1.0
           else if(hy_geometry == SPHERICAL) then
              lf = (xcent(i) - 0.5*dx)**2 / xcent(i)**2
              rf = (xcent(i) + 0.5*dx)**2 / xcent(i)**2
           else  !  CYLINDRICAL or POLAR
              lf = (xcent(i) - 0.5*dx) / xcent(i)
              rf = (xcent(i) + 0.5*dx) / xcent(i)
           end if

           if(hy_geometry == SPHERICAL .AND. NDIM > 1) then
              lfj = sin(ycent(j) - 0.5*dy) / sin(ycent(j))
              rfj = sin(ycent(j) + 0.5*dy) / sin(ycent(j))
           else
              lfj = 1.0
              rfj = 1.0
           end if



           ! Compute change in cell internal energy over the timestep:
           dutot =  utot_new - utot_old

           ! Compute change in ion/electron/radiation internal
           ! energies due to advection:
           durad_adv = - (dt/dx)*(rf*xflux(HY_ERAD_FLUX,i+1,j,    k    ) - lf*xflux(HY_ERAD_FLUX,i,j,k)) &
                       - (dt/dy)*(rfj*yflux(HY_ERAD_FLUX,i,  j+K2D,k    ) -lfj*yflux(HY_ERAD_FLUX,i,j,k)) &
                       - (dt/dz)*(   zflux(HY_ERAD_FLUX,i,  j,    k+K3D) -    zflux(HY_ERAD_FLUX,i,j,k))

           duele_adv = - (dt/dx)*(rf*xflux(HY_EELE_FLUX,i+1,j,    k    ) - lf*xflux(HY_EELE_FLUX,i,j,k)) &
                       - (dt/dy)*(rfj*yflux(HY_EELE_FLUX,i,  j+K2D,k    ) -lfj*yflux(HY_EELE_FLUX,i,j,k)) &
                       - (dt/dz)*(   zflux(HY_EELE_FLUX,i,  j,    k+K3D) -    zflux(HY_EELE_FLUX,i,j,k))
           
           duion_adv = - (dt/dx)*(rf*xflux(HY_EION_FLUX,i+1,j,    k    ) - lf*xflux(HY_EION_FLUX,i,j,k)) &
                       - (dt/dy)*(rfj*yflux(HY_EION_FLUX,i,  j+K2D,k    ) -lfj*yflux(HY_EION_FLUX,i,j,k)) &
                       - (dt/dz)*(   zflux(HY_EION_FLUX,i,  j,    k+K3D) -    zflux(HY_EION_FLUX,i,j,k))


           dueleAdvPlus = 0.0
           duionAdvPlus = 0.0
           duradAdvPlus = 0.0

#ifdef FLASH_USM_MHD
#ifdef FLASH_UHD_3T
           if (hy_hallVelocity) then
             dueleAdvPlus = - (dt/dx)*(  rf*( - Jp(1,i,j,    k    )) &
                                       - lf*( - Jm(1,i,j,    k    ))) &
                         - (dt/dy)*(     - rfj*Jp(2,i,j+K2D,k    ) &
                                         + lfj*Jm(2,i,j    ,k    )) &
                         - (dt/dz)*(     - Jp(3,i,j,    k+K3D) &
                                         + Jm(3,i,j,    k    ))
           endif
#endif
#endif

           lam3Rad = U(FLLM_VAR,i,j,k) * hy_lam3ScaleFactor
           lam3RadFloored = max(lam3Rad,hy_mtPresRatLambda3Min)

#define PREDICT_LAMBDA3_FOR_RADACCEL_FROM_PRESRAT
#define PREDICT_LAMBDA3_FOR_RADWORK_FROM_PRESRAT

           presFL=xflux(HY_P_FLUX   ,i,j,k); presFR=xflux(HY_P_FLUX   ,i+1,j,    k    )
           pmatFL=xflux(HY_PMAT_FLUX,i,j,k); pmatFR=xflux(HY_PMAT_FLUX,i+1,j,    k    )
           pradFL=xflux(HY_PRAD_FLUX,i,j,k); pradFR=xflux(HY_PRAD_FLUX,i+1,j,    k    )
           velxFL=lf*xflux(HY_VOLU_FLUX,i,j,k); velxFR=rf*xflux(HY_VOLU_FLUX,i+1,j,    k    )
           presGL=yflux(HY_P_FLUX   ,i,j,k); presGR=yflux(HY_P_FLUX   ,i  ,j+K2D,k    )
           pmatGL=yflux(HY_PMAT_FLUX,i,j,k); pmatGR=yflux(HY_PMAT_FLUX,i  ,j+K2D,k    )
           pradGL=yflux(HY_PRAD_FLUX,i,j,k); pradGR=yflux(HY_PRAD_FLUX,i  ,j+K2D,k    )
           velyGL=lfj*yflux(HY_VOLU_FLUX,i,j,k); velyGR=rfj*yflux(HY_VOLU_FLUX,i  ,j+K2D,k    )
           presHL=zflux(HY_P_FLUX   ,i,j,k); presHR=zflux(HY_P_FLUX   ,i  ,j   , k+K3D)
           pmatHL=zflux(HY_PMAT_FLUX,i,j,k); pmatHR=zflux(HY_PMAT_FLUX,i  ,j   , k+K3D)
           pradHL=zflux(HY_PRAD_FLUX,i,j,k); pradHR=zflux(HY_PRAD_FLUX,i  ,j   , k+K3D)
           velzHL=zflux(HY_VOLU_FLUX,i,j,k); velzHR=zflux(HY_VOLU_FLUX,i  ,j   , k+K3D)

#ifdef PREDICT_LAMBDA3_FOR_RADACCEL_FROM_PRESRAT
           if (scaleLorentz .NE. 0.0 .and. lam3Rad .NE. 0.0) then
#endif
              uDotGradPradPhys = &
                  (0.5/dx)*(xflux(HY_PRAD_FLUX,i+1,j,    k    )-xflux(HY_PRAD_FLUX,i,j,k)) &
                *(rf*xflux(HY_VOLU_FLUX,i+1,j,    k    )+lf*xflux(HY_VOLU_FLUX,i,j,k)) &
                + (0.5/dy)*(yflux(HY_PRAD_FLUX,i,  j+K2D,k    )-yflux(HY_PRAD_FLUX,i,j,k)) &
                *(rfj*yflux(HY_VOLU_FLUX,i,  j+K2D,k    )+lfj*yflux(HY_VOLU_FLUX,i,j,k)) &
                + (0.5/dz)*(zflux(HY_PRAD_FLUX,i,  j,    k+K3D)-zflux(HY_PRAD_FLUX,i,j,k)) &
                *(zflux(HY_VOLU_FLUX,i,  j,    k+K3D)+zflux(HY_VOLU_FLUX,i,j,k))
#ifdef PREDICT_LAMBDA3_FOR_RADACCEL_FROM_PRESRAT
           end if
#endif

#ifndef PREDICT_LAMBDA3_FOR_RADACCEL_FROM_PRESRAT
           duradAccel = lam3RadFloored * uDotGradPradPhys
#else
           lam3L = max(1.d-199,(presFL - pmatFL)) / max(pradFL,1.d-199)
           lam3R = max(1.d-199,(presFR - pmatFR)) / max(pradFR,1.d-199)
           lam3Lf = max(lam3L,hy_mtPresRatLambda3Min)
           lam3Rf = max(lam3R,hy_mtPresRatLambda3Min)
           lam3uDotGradPradPhys = &
                (0.5/dx)*(xflux(HY_PRAD_FLUX,i+1,j,    k    )-xflux(HY_PRAD_FLUX,i,j,k)) &
                *(velxFR*lam3Rf + velxFL*lam3Lf)
           if (NDIM > 1) then
              lam3L = max(1.d-199,(presGL - pmatGL)) / max(pradGL,1.d-199)
              lam3R = max(1.d-199,(presGR - pmatGR)) / max(pradGR,1.d-199)
              lam3Lf = max(lam3L,hy_mtPresRatLambda3Min)
              lam3Rf = max(lam3R,hy_mtPresRatLambda3Min)
              lam3uDotGradPradPhys = lam3uDotGradPradPhys &
                   + (0.5/dy)*(yflux(HY_PRAD_FLUX,i,  j+K2D,k    )-yflux(HY_PRAD_FLUX,i,j,k)) &
                    *(velyGR*lam3Rf + velyGL*lam3Lf)
              if (NDIM > 2) then
                 lam3L = max(1.d-199,(presHL - pmatHL)) / max(pradHL,1.d-199)
                 lam3R = max(1.d-199,(presHR - pmatHR)) / max(pradHR,1.d-199)
                 lam3Lf = max(lam3L,hy_mtPresRatLambda3Min)
                 lam3Rf = max(lam3R,hy_mtPresRatLambda3Min)
                 lam3uDotGradPradPhys = lam3uDotGradPradPhys &
                      + (0.5/dz)*(zflux(HY_PRAD_FLUX,i,  j,    k+K3D)-zflux(HY_PRAD_FLUX,i,j,k)) &
                       *(velzHR*lam3Rf + velzHL*lam3Lf)
              end if
           end if

           duradAccel = lam3uDotGradPradPhys
#endif

           duradAccel = dt * fB2r *       duradAccel

           uDotGradPmat = &
                  (0.5/dx)*(xflux(HY_PMAT_FLUX,i+1,j,    k    )-xflux(HY_PMAT_FLUX,i,j,k)) &
                *(rf*xflux(HY_VOLU_FLUX,i+1,j,    k    )+lf*xflux(HY_VOLU_FLUX,i,j,k)) &
                + (0.5/dy)*(yflux(HY_PMAT_FLUX,i,  j+K2D,k    )-yflux(HY_PMAT_FLUX,i,j,k)) &
                *(rfj*yflux(HY_VOLU_FLUX,i,  j+K2D,k    )+lfj*yflux(HY_VOLU_FLUX,i,j,k)) &
                + (0.5/dz)*(zflux(HY_PMAT_FLUX,i,  j,    k+K3D)-zflux(HY_PMAT_FLUX,i,j,k)) &
                *(zflux(HY_VOLU_FLUX,i,  j,    k+K3D)+zflux(HY_VOLU_FLUX,i,j,k))

           dueleAccel =           uDotGradPmat
           dueleAccel = dt * fB2  *       dueleAccel

           duionAccel = 0.0
           uion_new = 0.0       !DEV: revise later

           dueleWork = &
                  (0.5/dx)*(xflux(HY_PMAT_FLUX,i+1,j,    k    )+xflux(HY_PMAT_FLUX,i,j,k)) &
                *(rf*xflux(HY_VOLU_FLUX,i+1,j,    k    )-lf*xflux(HY_VOLU_FLUX,i,j,k)) &
                + (0.5/dy)*(yflux(HY_PMAT_FLUX,i,  j+K2D,k    )+yflux(HY_PMAT_FLUX,i,j,k)) &
                *(rfj*yflux(HY_VOLU_FLUX,i,  j+K2D,k    )-lfj*yflux(HY_VOLU_FLUX,i,j,k)) &
                + (0.5/dz)*(zflux(HY_PMAT_FLUX,i,  j,    k+K3D)+zflux(HY_PMAT_FLUX,i,j,k)) &
                *(zflux(HY_VOLU_FLUX,i,  j,    k+K3D)-zflux(HY_VOLU_FLUX,i,j,k))

           duionWork = 0.0
#ifndef PREDICT_LAMBDA3_FOR_RADWORK_FROM_PRESRAT
           duradWork = &
                  (0.5/dx)*(xflux(HY_PRAD_FLUX,i+1,j,    k    )+xflux(HY_PRAD_FLUX,i,j,k)) &
                *(rf*xflux(HY_VOLU_FLUX,i+1,j,    k    )-lf*xflux(HY_VOLU_FLUX,i,j,k)) &
                + (0.5/dy)*(yflux(HY_PRAD_FLUX,i,  j+K2D,k    )+yflux(HY_PRAD_FLUX,i,j,k)) &
                *(rfj*yflux(HY_VOLU_FLUX,i,  j+K2D,k    )-lfj*yflux(HY_VOLU_FLUX,i,j,k)) &
                + (0.5/dz)*(zflux(HY_PRAD_FLUX,i,  j,    k+K3D)+zflux(HY_PRAD_FLUX,i,j,k)) &
                *(zflux(HY_VOLU_FLUX,i,  j,    k+K3D)-zflux(HY_VOLU_FLUX,i,j,k))
           duradWork = duradWork * lam3RadFloored
#endif

           lam3L = max(1.d-199,(presFL - pmatFL)) / max(pradFL,1.d-199)
           lam3R = max(1.d-199,(presFR - pmatFR)) / max(pradFR,1.d-199)
           lam3Lf = max(lam3L,hy_mtPresRatLambda3Min)
           lam3Rf = max(lam3R,hy_mtPresRatLambda3Min)

#ifdef PREDICT_LAMBDA3_FOR_RADWORK_FROM_PRESRAT
           pradEffFL = presFL - pmatFL     ; pradEffFR = presFR - pmatFR
           pradEffGL = presGL - pmatGL     ; pradEffGR = presGR - pmatGR
           pradEffHL = presHL - pmatHL     ; pradEffHR = presHR - pmatHR
           pradEffFL = max(pradEffFL,pradFL*hy_mtPresRatLambda3Min)
           pradEffFR = max(pradEffFR,pradFR*hy_mtPresRatLambda3Min)
!!$           duradWork = &
!!$                (0.5/dx)*(lam3R*xflux(HY_PRAD_FLUX,i+1,j,k)+lam3L*xflux(HY_PRAD_FLUX,i,j,k)) &
!!$                *(rf*xflux(HY_VOLU_FLUX,i+1,j,k)-lf*xflux(HY_VOLU_FLUX,i,j,k))
           duradWork = &
                (0.5/dx)*(pradEffFR + pradEffFL) &
                *(rf*xflux(HY_VOLU_FLUX,i+1,j,k)-lf*xflux(HY_VOLU_FLUX,i,j,k))
#endif
           duradWorkPlus = 1.0/dx*  (fB1r+fB0r)*(lam3Rf-lam3Lf)*0.25*(velxFL+velxFR)&
                                     *(pradFL+pradFR)
           if (NDIM > 1) then
              lam3L = max(1.d-199,(presGL - pmatGL)) / max(pradGL,1.d-199)
              lam3R = max(1.d-199,(presGR - pmatGR)) / max(pradGR,1.d-199)
              lam3Lf = max(lam3L,hy_mtPresRatLambda3Min)
              lam3Rf = max(lam3R,hy_mtPresRatLambda3Min)
#ifdef PREDICT_LAMBDA3_FOR_RADWORK_FROM_PRESRAT
              pradEffGL = max(pradEffGL,pradGL*hy_mtPresRatLambda3Min)
              pradEffGR = max(pradEffGR,pradGR*hy_mtPresRatLambda3Min)
!!$              duradWork = duradWork  &
!!$                + (0.5/dy)*(lam3R*yflux(HY_PRAD_FLUX,i,j+K2D,k)+lam3L*yflux(HY_PRAD_FLUX,i,j,k)) &
!!$                 *(rfj*yflux(HY_VOLU_FLUX,i,j+K2D,k)-lfj*yflux(HY_VOLU_FLUX,i,j,k))
              duradWork = duradWork  &
                + (0.5/dy)*(pradEffGR + pradEffGL) &
                 *(rfj*yflux(HY_VOLU_FLUX,i,j+K2D,k)-lfj*yflux(HY_VOLU_FLUX,i,j,k))
#endif
              duradWorkPlus = duradWorkPlus + 1.0/dy*  (fB1r+fB0r)*(lam3Rf-lam3Lf)*0.25*(velyGL+velyGR)&
                                                        *(pradGL+pradGR)
              if (NDIM > 2) then
                 lam3L = max(1.d-199,(presHL - pmatHL)) / max(pradHL,1.d-199)
                 lam3R = max(1.d-199,(presHR - pmatHR)) / max(pradHR,1.d-199)
                 lam3Lf = max(lam3L,hy_mtPresRatLambda3Min)
                 lam3Rf = max(lam3R,hy_mtPresRatLambda3Min)
#ifdef PREDICT_LAMBDA3_FOR_RADWORK_FROM_PRESRAT
                 pradEffHL = max(pradEffHL,pradHL*hy_mtPresRatLambda3Min)
                 pradEffHR = max(pradEffHR,pradHR*hy_mtPresRatLambda3Min)
!!$                 duradWork = duradWork  &
!!$                      + (0.5/dz)*(lam3R*zflux(HY_PRAD_FLUX,i,j,k+K3D)+lam3L*zflux(HY_PRAD_FLUX,i,j,k)) &
!!$                       *(zflux(HY_VOLU_FLUX,i,j,k+K3D)-zflux(HY_VOLU_FLUX,i,j,k))
                 duradWork = duradWork  &
                      + (0.5/dz)*(pradEffHR + pradEffHL) &
                       *(zflux(HY_VOLU_FLUX,i,j,k+K3D)-zflux(HY_VOLU_FLUX,i,j,k))
#endif
                 duradWorkPlus = duradWorkPlus + 1.0/dz*  (fB1r+fB0r)*(lam3Rf-lam3Lf)*0.25*(velzHL+velzHR)&
                                                           *(pradHL+pradHR)
              end if
           end if

           duionWorkLike = duionAccel + duionWork
           dueleWorkLike = dueleAccel - dt * (fB1 +fB0)  *  dueleWork 
           duradWorkLike = duradAccel - dt * (fB1r+fB0r) * (duradWork + duradWorkPlus)

           ! for informative error messages:
           work_ele =  - dt       * dueleWork
           dueleWork = - dt * fB1 * dueleWork

           duradWorkPlus = - dt * duradWorkPlus
           work_rad =  - dt        * duradWork
           duradWork = - dt * fB1r * duradWork  +  (fB1r+fB0r) * duradWorkPlus

           ! Compute the total advected interal energy. There are two ways to do this:
           !   1) Directly from HY_EINT_FLUX
           !   2) As the sum of the electron, ion, radiation internal energies.
           !
           ! In a perfect world, these would be the same. However -
           ! they are not. Since we are scaling the component
           ! energies, I think it is more appropriate to go with
           ! option (2)
           !
           ! This is how you would use HY_EINT_FLUX directly if you want...
           ! dutot_adv = - (dt/dx)*(rf*xflux(HY_EINT_FLUX,i+1,j,    k    ) - lf*xflux(HY_EINT_FLUX,i,j,k)) &
           !             - (dt/dy)*(rfj*yflux(HY_EINT_FLUX,i,  j+K2D,k    ) -lfj*yflux(HY_EINT_FLUX,i,j,k)) &
           !             - (dt/dz)*(   zflux(HY_EINT_FLUX,i,  j,    k+K3D) -    zflux(HY_EINT_FLUX,i,j,k))

           dutot_adv = durad_adv + duion_adv + duele_adv
           dutotAdvPlus = dutot_adv                       + duradAdvPlus + duionAdvPlus + dueleAdvPlus

           ! Compute advected ion/electron/radiation internal energies:
           uele_adv =  uele_old + duele_adv + dueleAdvPlus
           uion_adv =  uion_old + duion_adv
           urad_adv =  urad_old + durad_adv

           duion = duion_adv                + scaleAccel*duionAccel + scaleWork*duionWork
           duele = duele_adv + dueleAdvPlus + scaleAccel*dueleAccel + scaleWork*dueleWork
           durad = durad_adv                + scaleAccel*duradAccel + scaleWork*duradWork


           if (hy_3TMode == HY3T_CASTROLIKE) then
              uele_new = uele_old + duele
              uion_new = uion_old + duion
              urad_new = urad_old + durad
              if( uele_new <= 0.0 .or. uion_new < 0.0 .or. urad_new < 0.0 ) then
                 call Logfile_stampMessage("[hy_uhd_unsplitUpdateCastroLike] Warning #1 in u{ele,ion,rad}_new:", .true.)
                 call adv_err_chk(.TRUE.)
                 call internal_shiftEints(duion, &
                      duele, &
                      durad, &
                      U(EION_VAR,i,j,k),U(EELE_VAR,i,j,k), &
                      U(ERAD_VAR,i,j,k), &
                      dens_old, U(DENS_VAR,i,j,k),densNewInv,i,j,k)
                 uele_new = uele_old + duele
                 uion_new = uion_old + duion
                 urad_new = urad_old + durad
                 call Logfile_stampMessage("[hy_uhd_unsplitUpdateCastroLike] Info #1a in u{ele,ion,rad}_new:", .true.)
                 call adv_err_chk(.TRUE.)
              end if
           end if

!!$           ! Check for negative advected internal energies:
!!$           if( uele_adv <= 0.0 .or. uion_adv < 0.0 .or. urad_adv < 0.0 ) then
!!$              call Logfile_stampMessage("[hy_uhd_unsplitUpdateCastroLike] Error in u{ele,ion,rad}_adv:", .true.)
!!$              call adv_err_chk(.TRUE.)
!!$              call internal_shiftEints(duion_adv, &
!!$                                       duele_adv, &
!!$                                       durad_adv, &
!!$                                       U(EION_VAR,i,j,k),U(EELE_VAR,i,j,k), &
!!$                                       U(ERAD_VAR,i,j,k), &
!!$                                       dens_old, U(DENS_VAR,i,j,k),densNewInv,i,j,k)
!!$              uele_adv =  uele_old + duele_adv
!!$              uion_adv =  uion_old + duion_adv
!!$              urad_adv =  urad_old + durad_adv
!!$              if( uele_adv <= 0.0 .or. uion_adv < 0.0 .or. urad_adv < 0.0 ) then
!!$                 call Logfile_stampMessage("[hy_uhd_unsplitUpdateCastroLike] Error REMAINS in u{ele,ion,rad}_adv:", .true.)
!!$                 call adv_err_chk(.TRUE.)
!!$              end if
!!$           end if

           !! Subtract Qohm from Eint and add it after ragelike to eele
           Qohm = 0.0              
#ifdef FLASH_USM_MHD
           if (hy_useMagneticResistivity) then
             Qohm = scrch_Ptr(HY_XN06_SCRATCHCTR_VAR,i,j,k)
           endif   
#endif

           if (hy_3TMode == HY3T_RAGELIKE) then
              call hy_uhd_ragelike(dutot-dutotAdvPlus-Qohm, U(:,i,j,k), dens_old, &
                   duele_adv+dueleAdvPlus, duion_adv, durad_adv, &
                   xcent(i), ycent(j), zcent(k), &
                   uele_new, uion_new, urad_new)
           else if (hy_3TMode == HY3T_CASTROLIKE) then
              inShock  = ( U(SHOK_VAR,i,j,k) .GT. 0.0 )
              if (do2T) then
                 uref = max(uele_old,uele_adv,uele_new,abs(dueleAccel),abs(dueleWork),abs(duele_adv))
              else
                 uref = max(uion_old,uion_adv,uion_new,abs(duionAccel),abs(duionWork),abs(duion_adv))
              end if
              if (inShock .AND. &
                  dutot .GE. (duele+duion+durad)+Qohm .AND. &
                  dutot - (duele+duion+durad+Qohm) .LE. &
                    1.0e5 * uref ) then
                 ! Add it all to the matter!
#if EOSCOMP_MATTER > 0
                 duDirect = min(5.0*uref, dutot - (duele+duion+durad+Qohm) )
                 duele = duele + duDirect
#  if (0)
                 uele_new = uele_old + duele
                 duRad = duRad + dutot - (duele+duion+durad+Qohm)
!!$                 !urad_new = urad_old + durad
                 urad_new = max(urad_old + durad , utot_new - (uele_new + uion_new + Qohm))
#  else
                 call localRageLike(dutot-(duele+duion+durad)-Qohm, U(:,i,j,k), dens_old, &
                   uele_old,  uion_old,  urad_old, &
                   duele,     duion,     durad , &
                   xcent(i), ycent(j), zcent(k), &
                   uele_new, uion_new, urad_new)
#  endif
#else
                 duDirect = min(1.0*uref, dutot - (duele+duion+durad+Qohm) ) * 0.999455
                 duion = duion + duDirect
                 call hy_uhd_ragelike(dutot-(duele+duion+durad)-Qohm, U(:,i,j,k), dens_old, &
                   duele,     duion,     durad , &
                   xcent(i), ycent(j), zcent(k), &
                   uele_new, uion_new, urad_new)
#endif
              else
                 call localRageLike(dutot-(duele+duion+durad)-Qohm, U(:,i,j,k), dens_old, &
                   uele_old,  uion_old,  urad_old, &
                   duele,     duion,     durad , &
                   xcent(i), ycent(j), zcent(k), &
                   uele_new, uion_new, urad_new, &
                   dueleWorkLike,duradWorkLike,duionWorkLike)
              end if
              if( uele_new <= 0.0 .or. uion_new < 0.0 .or. urad_new < 0.0 ) then
                 call Logfile_stampMessage("[hy_uhd_unsplitUpdateCastroLike] Warning #8 in u{ele,ion,rad}_new:", .true.)
                 call adv_err_chk(.TRUE.)
              end if
              call internal_shiftEints(uion_new, &
                                       uele_new, &
                                       urad_new, &
                                       0.0, 0.0, 0.0, &
                                       1.0, U(DENS_VAR,i,j,k),densNewInv,i,j,k)
           else
              uele_new = uele_adv
              uion_new = uion_adv
              urad_new = urad_adv
           end if

           if (scaleLorentz .NE. 0.0 .and. lam3Rad .NE. 0.0 .and. uDotGradPradPhys .NE. 0.0) then
              call Opacity(U(:,i,j,k),1,kappaP,kappaDummy,kappaR)
              kappaRatio = kappaP / kappaR

              if (kappaRatio .NE. 0.0) then
                 duradLamPlus = 2.0*kappaRatio * lam3Rad * uDotGradPradPhys
                 scaledDuradLamPlus = -dt * scaleLorentz * duradLamPlus

                 urad_new = urad_new + scaledDuradLamPlus
                 uele_new = uele_new - scaledDuradLamPlus
                 if( uele_new <= 0.0 .or. uion_new < 0.0 .or. urad_new < 0.0 ) then
                    call Logfile_stampMessage("[hy_uhd_unsplitUpdateCastroLike] Warning #9 in u{ele,ion,rad}_new:", .true.)
                    call adv_err_chk(.TRUE.)
                 end if
                 call internal_shiftEints(uion_new, &
                                          uele_new, &
                                          urad_new, &
                                          0.0, 0.0, 0.0, &
                                          1.0, U(DENS_VAR,i,j,k),densNewInv,i,j,k)
              end if
           end if


           !! Update Uele with Qohm (it is actually dens*Qohm*dt at t = n)
           uele_new = uele_new + Qohm
           
           U(EELE_VAR,i,j,k) = uele_new * densNewInv
           U(EION_VAR,i,j,k) = uion_new * densNewInv
           U(ERAD_VAR,i,j,k) = urad_new * densNewInv           

           ! The values of eion, eele, or erad may be negative or zero
           ! at this point for unknown reasons. Check for this:
           if( U(EELE_VAR, i,j,k) <= 0.0 .or. &
               U(EION_VAR, i,j,k) <  0.0 .or. &
               U(ERAD_VAR, i,j,k) <  0.0 .or. &
               uele_new <= 0.0 .or. &
               uion_new <  0.0 .or. &
               urad_new <  0.0 ) then
               
              call Logfile_stampMessage("[hy_uhd_unsplitUpdateCastroLike] ERROR DETECTED", .true.)
              call Logfile_stampMessage("  Zero or negative advected specific internal energy detected", .true.)

              write(errmsg, "(a,1pe20.13)") "  DENS_OLD = ", dens_old
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "  DENS_NEW = ", dens_new
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "  DUTOT_ADV = ", dutot_adv
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "  DUTOT_NEW = ", dutot
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "  UTOT_OLD = ", utot_old
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "  UTOT_NEW = ", utot_new
              call Logfile_stampMessage(trim(errmsg), .true.)

              call Logfile_stampMessage("  OLD INTERNAL ENERGY DENSITIES:")
              write(errmsg, "(a,1pe20.13)") "    UELE_OLD = ", uele_old
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "    UION_OLD = ", uion_old
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "    URAD_OLD = ", urad_old
              call Logfile_stampMessage(trim(errmsg), .true.)

              call Logfile_stampMessage("  NEW INTERNAL ENERGY DENSITIES:")              
              write(errmsg, "(a,1pe20.13)") "    UELE_NEW = ", U(EELE_VAR,i,j,k) * U(DENS_VAR,i,j,k)
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "    UION_NEW = ", U(EION_VAR,i,j,k) * U(DENS_VAR,i,j,k)
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "    URAD_NEW = ", U(ERAD_VAR,i,j,k) * U(DENS_VAR,i,j,k)
              call Logfile_stampMessage(trim(errmsg), .true.)
              
              write(errmsg, "(a,1pe20.13)") "    duele     = ", duele
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "    duion     = ", duion
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "    durad     = ", durad
              call Logfile_stampMessage(trim(errmsg), .true.)

              call Logfile_stampMessage("  NEW SPECIFIC INTERNAL ENERGIES:")
              write(errmsg, "(a,1pe20.13)") "    EELE_VAR = ", U(EELE_VAR,i,j,k)
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "    EION_VAR = ", U(EION_VAR,i,j,k)
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "    ERAD_VAR = ", U(ERAD_VAR,i,j,k)
              call Logfile_stampMessage(trim(errmsg), .true.)

              call Logfile_stampMessage("  ADVECTED INTERNAL ENERGY DENSITIES:")
              write(errmsg, "(a,1pe20.13)") "    uele_adv  = ", uele_adv
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "    uion_adv  = ", uion_adv
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "    urad_adv  = ", urad_adv
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "    duele_adv = ", duele_adv
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "    duion_adv = ", duion_adv
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "    durad_adv = ", durad_adv
              call Logfile_stampMessage(trim(errmsg), .true.)
              
              write(errmsg, "(a,1pe20.13)") "  PELE_VAR = ", U(PELE_VAR,i,j,k)
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "  PION_VAR = ", U(PION_VAR,i,j,k)
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "  PRAD_VAR = ", U(PRAD_VAR,i,j,k)
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "  PELE_ADV = ", pele_adv
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "  PION_ADV = ", pion_adv
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "  PRAD_ADV = ", prad_adv
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "  TELE_VAR = ", U(TELE_VAR,i,j,k)
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "  TION_VAR = ", U(TION_VAR,i,j,k)
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "  TRAD_VAR = ", U(TRAD_VAR,i,j,k)
              call Logfile_stampMessage(trim(errmsg), .true.)
              
              write(errmsg, "(a,1pe20.13)") "  CELL X = ", xcent(i)
              call Logfile_stampMessage(trim(errmsg), .true.)
              
              write(errmsg, "(a,1pe20.13)") "  CELL Y = ", ycent(j)
              call Logfile_stampMessage(trim(errmsg), .true.)
              
              write(errmsg, "(a,1pe20.13)") "  CELL Z = ", zcent(k)
              call Logfile_stampMessage(trim(errmsg), .true.)

              call Logfile_stampMessage("  You might be able to get things to keep running by trying", .true.)
              call Logfile_stampMessage("  more diffusive hydro options (reduce CFL, order, etc...)", .true.)

              if(uion_adv < 0.0) then
                 call Logfile_stampMessage("", .true.)
                 call Logfile_stampMessage("NEGATIVE ADVECTED ION ENERGY:", .true.)
                 call flux_error(HY_EION_FLUX)
              end if

              if(uele_adv < 0.0) then
                 call Logfile_stampMessage("", .true.)
                 call Logfile_stampMessage("NEGATIVE ADVECTED ELECTRON ENERGY:", .true.)
                 call flux_error(HY_EELE_FLUX)
              end if

              if(urad_adv < 0.0) then
                 call Logfile_stampMessage("", .true.)
                 call Logfile_stampMessage("NEGATIVE ADVECTED RADIATION ENERGY:", .true.)
                 call flux_error(HY_ERAD_FLUX)
              end if

              call Logfile_stampMessage("  WORK TERMS:")
              call computeAddtlDiagTerms()

              write(errmsg, "(a,1pe20.13)") "    divv = ", divv
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "    pele = ", pele
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "    prad = ", prad
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "    work_ele = ", work_ele
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "    work_rad = ", work_rad
              call Logfile_stampMessage(trim(errmsg), .true.)

              call Driver_abortFlash("[hy_uhd_unsplitUpdateCastroLike] Negative 3T internal energy, CHECK LOG")
              
           end if
           
        end do
     end do
  end do

  ! Release block pointers
  call Grid_releaseBlkPtr(blockID,U,CENTER)

  ! Deallocate arrays:
  deallocate(xcent)
  deallocate(ycent)
  deallocate(zcent)

contains

  subroutine localRageLike(uextra, soln, dens_old, &
     uele_old, uion_old, urad_old, &
     duele_adv, duion_adv, durad_adv, &
     xc, yc, zc, &
     uele_new, uion_new, urad_new, &
     extraEle, extraRad, extraIon)
    use hy_uhd_MultiTempData, ONLY: hy_3Ttry_D, hy_3Ttry_G

  ! Arguments:
    real, intent(in)  :: uextra
    real, intent(in)  :: soln(NUNK_VARS)
    real, intent(in)  :: dens_old
    real, intent(in)  :: uele_old
    real, intent(in)  :: uion_old
    real, intent(in)  :: urad_old
    real, intent(in)  :: duele_adv
    real, intent(in)  :: duion_adv
    real, intent(in)  :: durad_adv
    real, intent(in)  :: xc 
    real, intent(in)  :: yc 
    real, intent(in)  :: zc
    real, intent(out) :: uele_new
    real, intent(out) :: uion_new
    real, intent(out) :: urad_new
    real, intent(in),OPTIONAL :: extraEle, extraRad, extraIon

    real :: avgPradPhys, avgPele, avgPion
    real :: uextraEle, uextraIon, uextraRad, uextraSum
    real :: uextraTmp, uIonTmp, uEleTmp, uRadTmp
    real :: duele, duion, durad

    if (uele_old+duele_adv   >  0.0 .AND. &
        uion_old+duion_adv .GE. 0.0 .AND. &
        urad_old+durad_adv .GE. 0.0 .AND. &
        hy_3Ttry_G .NE. 1) then
       call hy_uhd_ragelike(uextra, soln, dens_old, &
            duele_adv, duion_adv, durad_adv, &
            xc, yc, zc, &
            uele_new, uion_new, urad_new, &
            stepFactor=0.25*real(hy_3Ttry_G))
       avgPradPhys = &
                  (xflux(HY_PRAD_FLUX,i+1,j,    k    )+xflux(HY_PRAD_FLUX,i,j,k))
       if (NDIM>1) then
          avgPradPhys = avgPradPhys  &
                + (yflux(HY_PRAD_FLUX,i,  j+K2D,k    )+yflux(HY_PRAD_FLUX,i,j,k))
          if (NDIM>2) then
             avgPradPhys = avgPradPhys  &
                  + (zflux(HY_PRAD_FLUX,i,  j,    k+K3D)+zflux(HY_PRAD_FLUX,i,j,k))
             avgPradPhys = avgPradPhys  / 6.0
          else
             avgPradPhys = avgPradPhys  * 0.25
          end if
       else
          avgPradPhys = avgPradPhys  * 0.5
       end if
       avgPele = &
                  (xflux(HY_PMAT_FLUX,i+1,j,    k    )+xflux(HY_PMAT_FLUX,i,j,k))
       if (NDIM>1) then
          avgPele = avgPele  &
                + (yflux(HY_PMAT_FLUX,i,  j+K2D,k    )+yflux(HY_PMAT_FLUX,i,j,k))
          if (NDIM>2) then
             avgPele = avgPele  &
                  + (zflux(HY_PMAT_FLUX,i,  j,    k+K3D)+zflux(HY_PMAT_FLUX,i,j,k))
             avgPele = avgPele  / 6.0
          else
             avgPele = avgPele  * 0.25
          end if
       else
          avgPele = avgPele  * 0.5
       end if
       avgPion = -1.0            !diag only
    else
       avgPradPhys = &
                  (xflux(HY_PRAD_FLUX,i+1,j,    k    )+xflux(HY_PRAD_FLUX,i,j,k))
       if (NDIM>1) then
          avgPradPhys = avgPradPhys  &
                + (yflux(HY_PRAD_FLUX,i,  j+K2D,k    )+yflux(HY_PRAD_FLUX,i,j,k))
          if (NDIM>2) then
             avgPradPhys = avgPradPhys  &
                  + (zflux(HY_PRAD_FLUX,i,  j,    k+K3D)+zflux(HY_PRAD_FLUX,i,j,k))
             avgPradPhys = avgPradPhys  / 6.0
          else
             avgPradPhys = avgPradPhys  * 0.25
          end if
       else
          avgPradPhys = avgPradPhys  * 0.5
       end if
       avgPele = &
                  (xflux(HY_PMAT_FLUX,i+1,j,    k    )+xflux(HY_PMAT_FLUX,i,j,k))
       if (NDIM>1) then
          avgPele = avgPele  &
                + (yflux(HY_PMAT_FLUX,i,  j+K2D,k    )+yflux(HY_PMAT_FLUX,i,j,k))
          if (NDIM>2) then
             avgPele = avgPele  &
                  + (zflux(HY_PMAT_FLUX,i,  j,    k+K3D)+zflux(HY_PMAT_FLUX,i,j,k))
             avgPele = avgPele  / 6.0
          else
             avgPele = avgPele  * 0.25
          end if
       else
          avgPele = avgPele  * 0.5
       end if
       avgPion = 0.0            !DEV: review later

       if (do2T) then
          avgPele = avgPele + avgPion
          avgPion = 0.0
       end if

       if (present(extraRad) .AND. present(extraEle)) then
          uextraRad = extraRad
          uextraEle = extraEle
          if (present(extraIon)) then
             uextraIon = extraIon
          else
             uextraIon = 0.0
          end if
          if (do2T) then
             uextraEle = uextraEle + uextraIon
             uextraIon = 0.0
          end if

          uextraTmp = uextra
          uextraSum = uextraIon + uextraEle + uextraRad
          uIonTmp = 0; uEleTmp = 0; uRadTmp = 0
          if (abs(uextraIon + uextraEle + uextraRad) < &
              abs(uextraIon) + abs(uextraEle) + abs(uextraRad) ) then

!!$             if (abs(uextra) > abs( (uextraIon + uextraEle)*(fB0 ) &
!!$                                   + uextraRad             *(fB0r))) then
                if (uextra*uextraSum > 0 .AND. &
                    abs(uextra) > abs(uextraSum)) then
                   if (uextra*uextraRad < 0.0) then
                         uextraTmp = uextra-uextraRad
                         uRadTmp = uextraRad; uextraRad = 0
                   end if
                   if (uextra*uextraIon < 0.0) then
                      uextraTmp = uextra-uextraIon
                      uIonTmp = uextraIon; uextraIon = 0
                   end if
                   if (uextra*uextraEle < 0.0) then
                      uextraTmp = uextra-uextraEle
                      uEleTmp = uextraEle; uextraEle = 0
                   end if
                else if (uextra*uextraSum < 0.0) then
                   if (uextra*uextraRad < 0.0) then
                      uextraRad = -uextraRad
                   end if
                   if (uextra*uextraIon < 0.0) then
                      uextraIon = -uextraIon
                   end if
                   if (uextra*uextraEle < 0.0) then
                      uextraEle = -uextraEle
                   end if
!!$                else
!!$                   if (uextra*uextraIon < 0.0) uextraIon = 0.0
!!$                   if (uextra*uextraEle < 0.0) uextraEle = 0.0
!!$                   if (uextra*uextraRad < 0.0) uextraRad = 0.0
!!$                   uextraRad = abs(extraRad)
!!$                   uextraEle = abs(extraEle)
!!$                   uextraIon = abs(extraIon)
                end if
!!$             end if

          else
             if (abs(uextraIon + uextraEle + uextraRad) .LE. TINY(uextraIon)) then
                uextraIon = 1.0 
             end if
          end if

#if(0)
          uextraRad = abs(extraRad)
          uextraEle = abs(extraEle)
          if (present(extraIon)) then
             uextraIon = abs(extraIon)
          else
             uextraIon = 0.0
          end if
#endif

          if (do2T) then
             uextraEle = uextraEle + uextraIon
             uextraIon = 0.0
          end if
          call Hydro_recalibrateEints(uextraTmp,uextraIon,uextraEle,uextraRad)
          uextraRad = uextraRad + uRadTmp
          uextraIon = uextraIon + uIonTmp
          uextraEle = uextraEle + uEleTmp
       else
          uextraIon = max(0.0,avgPion)
          uextraEle = max(0.0,avgPele)
          uextraRad = max(0.0,avgPradPhys)
          if (uextraIon + uextraEle + uextraRad .LE. 0.0) uextraIon = 1.0 

          if (do2T) then
             uextraEle = uextraEle + uextraIon
             uextraIon = 0.0
          end if

          call Hydro_recalibrateEints(uextra,uextraIon,uextraEle,uextraRad)
       end if


       duele = duele_adv + uextraEle
       duion = duion_adv + uextraIon
       durad = durad_adv + uextraRad

       uele_new = uele_old + (duele_adv + uextraEle)
       uion_new = uion_old + (duion_adv + uextraIon)
       urad_new = urad_old + (durad_adv + uextraRad)

!!$       call internal_shiftEints(duion, &
!!$            duele, &
!!$            durad, &
!!$            soln(EION_VAR),soln(EELE_VAR), &
!!$            soln(ERAD_VAR), &
!!$            dens_old, soln(DENS_VAR),densNewInv,i,j,k)
       uele_new = uele_old + duele
       uion_new = uion_old + duion
       urad_new = urad_old + durad

    end if

    pele_adv = avgPele
    pion_adv = avgPion
    prad_adv = avgPradPhys

  end subroutine localRageLike

  subroutine computeAddtlDiagTerms()
    divv = &
                  (1.0/dx)* &
                 (rf*xflux(HY_VOLU_FLUX,i+1,j,    k    )-lf*xflux(HY_VOLU_FLUX,i,j,k)) &
                + (1.0/dy)* &
                 (rfj*yflux(HY_VOLU_FLUX,i,  j+K2D,k    )-lfj*yflux(HY_VOLU_FLUX,i,j,k)) &
                + (1.0/dz)* &
                 (zflux(HY_VOLU_FLUX,i,  j,    k+K3D)-zflux(HY_VOLU_FLUX,i,j,k))

    work_rad = &
                  (0.5/dx)*(xflux(HY_PRAD_FLUX,i+1,j,    k    )+xflux(HY_PRAD_FLUX,i,j,k)) &
                *(rf*xflux(HY_VOLU_FLUX,i+1,j,    k    )-lf*xflux(HY_VOLU_FLUX,i,j,k)) &
                + (0.5/dy)*(yflux(HY_PRAD_FLUX,i,  j+K2D,k    )+yflux(HY_PRAD_FLUX,i,j,k)) &
                *(rfj*yflux(HY_VOLU_FLUX,i,  j+K2D,k    )-lfj*yflux(HY_VOLU_FLUX,i,j,k)) &
                + (0.5/dz)*(zflux(HY_PRAD_FLUX,i,  j,    k+K3D)+zflux(HY_PRAD_FLUX,i,j,k)) &
                *(zflux(HY_VOLU_FLUX,i,  j,    k+K3D)-zflux(HY_VOLU_FLUX,i,j,k))
    work_rad = - dt * work_rad

    if (uele_old+duele   >  0.0 .AND. &
        uion_old+duion .GE. 0.0 .AND. &
        urad_old+durad .GE. 0.0) then
       call hy_uhd_getPressure( &
            (uele_old+duele) / dens_new, &
            (uion_old+duion) / dens_new, &
            (urad_old+durad) / dens_new, &
            U(:,i,j,k), &
            pele, &
            pion_adv, &
            prad)
    end if

  end subroutine computeAddtlDiagTerms


  subroutine adv_err_chk(cont)
    implicit none
    logical,optional,intent(in) :: cont

    call Logfile_stampMessage("[hy_uhd_unsplitUpdateCastroLike] ERROR DETECTED (adv_err_chk)", .true.)
    call Logfile_stampMessage("  Zero or negative specific internal energy detected", .true.)

    write(errmsg, "(a,1pe20.13)") "  DENS_OLD = ", dens_old
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "  DENS_NEW = ", dens_new
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "  DUTOT_ADV = ", dutot_adv
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "  DUTOT_NEW = ", dutot
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "  UTOT_OLD = ", utot_old
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "  UTOT_NEW = ", utot_new
    call Logfile_stampMessage(trim(errmsg), .true.)

    call Logfile_stampMessage("  OLD INTERNAL ENERGY DENSITIES:")
    write(errmsg, "(a,1pe20.13)") "    UELE_OLD = ", uele_old
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "    UION_OLD = ", uion_old
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "    URAD_OLD = ", urad_old
    call Logfile_stampMessage(trim(errmsg), .true.)

    call Logfile_stampMessage("  NEW INTERNAL ENERGY DENSITIES:")
    write(errmsg, "(a,1pe20.13)") "    UELE_NEW = ", uele_new
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "    UION_NEW = ", uion_new
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "    URAD_NEW = ", urad_new
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "    added duele= ", uele_new-uele_old-duele
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "    added durad= ", urad_new-urad_old-durad
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "    duele     = ", duele
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "    duion     = ", duion
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "    durad     = ", durad
    call Logfile_stampMessage(trim(errmsg), .true.)

    if (hy_3Ttry_G ==1) then
       write(errmsg, "(a,1pe20.13)") "    dueleWkLike= ", dueleWorkLike
       call Logfile_stampMessage(trim(errmsg), .true.)

       write(errmsg, "(a,1pe20.13)") "    duionWkLike= ", duionWorkLike
       call Logfile_stampMessage(trim(errmsg), .true.)

       write(errmsg, "(a,1pe20.13)") "    duradWkLike= ", duradWorkLike
       call Logfile_stampMessage(trim(errmsg), .true.)
    end if

    call Logfile_stampMessage("  ADVECTED INTERNAL ENERGY DENSITIES:")
    write(errmsg, "(a,1pe20.13)") "    uele_adv  = ", uele_adv
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "    uion_adv  = ", uion_adv
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "    urad_adv  = ", urad_adv
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "    duele_adv = ", duele_adv
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "    duion_adv = ", duion_adv
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "    durad_adv = ", durad_adv
    call Logfile_stampMessage(trim(errmsg), .true.)

    if (hy_3Ttry_G ==1) then
              write(errmsg, "(a,1pe20.13)") "  PELE_VAR = ", U(PELE_VAR,i,j,k)
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "  PION_VAR = ", U(PION_VAR,i,j,k)
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "  PRAD_VAR = ", U(PRAD_VAR,i,j,k)
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "  PELE_ADV = ", pele_adv
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "  PION_ADV = ", pion_adv
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "  PRAD_ADV = ", prad_adv
              call Logfile_stampMessage(trim(errmsg), .true.)

              write(errmsg, "(a,1pe20.13)") "  lam3Rad  = ", lam3Rad
              call Logfile_stampMessage(trim(errmsg), .true.)
    end if

    call Logfile_stampMessage("  WORK TERMS:")
    call computeAddtlDiagTerms()

    write(errmsg, "(a,1pe20.13)") "    divv = ", divv
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "    pele = ", pele
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "    prad = ", prad
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "    dueleWork = ", dueleWork
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "    duradWork = ", duradWork
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "    dueleAccel = ", dueleAccel
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "    duradAccel = ", duradAccel
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "    work_ele = ", work_ele
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "    work_rad = ", work_rad
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "    duradWork+ =", duradWorkPlus
    call Logfile_stampMessage(trim(errmsg), .true.)


    write(errmsg, "(a,1pe20.13)") "  TELE_VAR = ", U(TELE_VAR,i,j,k)
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "  TION_VAR = ", U(TION_VAR,i,j,k)
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "  TRAD_VAR = ", U(TRAD_VAR,i,j,k)
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "  CELL X = ", xcent(i)
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "  CELL Y = ", ycent(j)
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1pe20.13)") "  CELL Z = ", zcent(k)
    call Logfile_stampMessage(trim(errmsg), .true.)

    call Logfile_stampMessage("  You might be able to get things to keep running by trying", .true.)
    call Logfile_stampMessage("  more diffusive hydro options (reduce CFL, order, etc...)", .true.)

    if(uion_adv < 0.0) then
       call Logfile_stampMessage("", .true.)
       call Logfile_stampMessage("NEGATIVE ADVECTED ION ENERGY:", .true.)
       call flux_error(HY_EION_FLUX)
    end if

    if(uele_adv < 0.0) then
       call Logfile_stampMessage("", .true.)
       call Logfile_stampMessage("NEGATIVE ADVECTED ELECTRON ENERGY:", .true.)
       call flux_error(HY_EELE_FLUX)
    end if

    if(urad_adv < 0.0) then
       call Logfile_stampMessage("", .true.)
       call Logfile_stampMessage("NEGATIVE ADVECTED RADIATION ENERGY:", .true.)
       call flux_error(HY_ERAD_FLUX)
    end if

    if (.NOT. cont) &
         call Driver_abortFlash("[hy_uhd_unsplitUpdateCastroLike] Negative 3T internal energy, CHECK LOG")
    
  end subroutine adv_err_chk


  subroutine flux_error(flux)
    implicit none
    integer, intent(in) :: flux

    write(errmsg, "(a,1p,2(e20.13):2G10.3)") "  Geometric Factors:", lf, rf, lfj, rfj
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,2(1pe20.13))") "  XFLUXES: ", xflux(flux,i,j,k), xflux(flux,i+1,j,k)
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,2(1pe20.13))") "  YFLUXES: ", yflux(flux,i,j,k), yflux(flux,i,j+K2D,k)
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,2(1pe20.13))") "  ZFLUXES: ", zflux(flux,i,j,k), zflux(flux,i,j,k+K3D)
    call Logfile_stampMessage(trim(errmsg), .true.)

    call Logfile_stampMessage("  CHANGE IN ENERGY DENSITY DUE TO FLUXES IN EACH DIRECTION:", .true.)

    write(errmsg, "(a,1(1pe20.13))") "  X DIR: ", &
         - (dt/dx)*(rf*xflux(flux,i+1,j,k) - lf*xflux(flux,i,j,k))
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1(1pe20.13))") "  Y DIR: ", &
         - (dt/dy)*(yflux(flux,i,j+K2D,k) - yflux(flux,i,j,k))
    call Logfile_stampMessage(trim(errmsg), .true.)

    write(errmsg, "(a,1(1pe20.13))") "  Z DIR: ", &
         - (dt/dz)*(zflux(flux,i,j,k+K3D) - zflux(flux,i,j,k))
    call Logfile_stampMessage(trim(errmsg), .true.)    

  end subroutine flux_error

  subroutine internal_shiftEints(deltaEion,deltaEele, &
                               deltaErad, &
                               eintIon_o,eintEle_o, &
                               eintRad_o, &
                               rho_o, newDens, inv_new_dens, i,j,k)

    use hy_uhd_MultiTempData, ONLY : hy_smallEion,hy_smallEele,hy_smallErad

    implicit none
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
    if (eAllNewAbove .LT. -0.5*spacing(hy_smallEion+hy_smallEele+hy_smallErad)) then
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
!#define DEBUG_MAR2015
#ifdef DEBUG_MAR2015
    if (oldRho .LE. 0.0 .OR. newRho .LE. 0.0 .OR. newRhoInv .LE. 0.0) &
         print*,'!!oldRho,newRho,newRhoInv:',oldRho,newRho,newRhoInv,'!!!!!'
!!$    print*,'e{All,Ion,Ele,Rad}NewAbove =',eAllNewAbove, eIonNewAbove, eEleNewAbove, eRadNewAbove,i,inShock
800 format(a,1P,3(G26.19),' @ (',i3,',',i3,',',i3,')')
    print 800,'_shift deltas BEFORE:',deltaEion,deltaEele,deltaErad,i,j,k
#endif

    eIonUnder = min(eintIon_o * oldRho + deltaEion - hy_smallEion*newRho, 0.0)
    eEleUnder = min(eintEle_o * oldRho + deltaEele - hy_smallEele*newRho, 0.0)
    eRadUnder = min(eintRad_o * oldRho + deltaErad - hy_smallErad*newRho, 0.0)

#ifdef DEBUG_MAR2015
    print 800,'_shift Unders BEFORE:',eIonUnder,eEleUnder,eRadUnder,i,j,k
#endif
    if (eIonUnder < 0.0) then
       deltaEion = deltaEion - eIonUnder
       if (eeleUnder == 0.0) then
#ifdef DEBUG_MAR2015
801       format('Shifting ',1P,G25.18,' from ',a,' to ',a,'...')
          print 801, -eIonUnder,'Electrons','Ions'
#endif
!!$           print '(a,2(G30.20))','deltaEele bef.:', deltaEele,eIonUnder
          deltaEele = deltaEele + eIonUnder !shift ion energy deficit to electrons
!!$           print '(a,2(G30.20))','deltaEele aft.:', deltaEele,spacing(deltaEele)
          eEleUnder = min(eintEle_o * oldRho + deltaEele - hy_smallEele*newRho, 0.0)
       else
#ifdef DEBUG_MAR2015
          print 801, -eIonUnder,'Radiation','Ions'
#endif
          deltaErad = deltaErad + eIonUnder !shift ion energy deficit to radiation
          eRadUnder = min(eintRad_o * oldRho + deltaErad - hy_smallErad*newRho, 0.0)
       end if
       eIonUnder = 0.0
    end if
    if (eEleUnder < 0.0) then  !last ditch effort:
#ifdef DEBUG_MAR2015
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
#ifdef DEBUG_MAR2015
          print 801, -eRadUnder,'Electrons','Radiation'
#endif
          deltaEele = deltaEele + eRadUnder !shift radiations energy deficit to electrons
!!$           eEleUnder = min(eintEle_o * oldRho + deltaEele - hy_smallEele*newRho, 0.0)
       else
#ifdef DEBUG_MAR2015
          print 801, -eRadUnder,'Ions','Radiation'
#endif
          deltaEion = deltaEion + eRadUnder !shift radiation energy deficit to ions
!!$           eIonUnder = min(eintIon_o * oldRho + deltaEion - hy_smallEion*newRho, 0.0)
       end if
    end if

#ifdef DEBUG_MAR2015
    print 800,'_shift deltas AFTER: ',deltaEion,deltaEele,deltaErad,i,j,k
#endif
  end subroutine internal_shiftEints

end Subroutine hy_uhd_unsplitUpdateCastroLike
