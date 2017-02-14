!!****if* source/physics/Hydro/HydroMain/unsplit/multiTemp/hy_uhd_unsplitUpdateMultiTemp
!!
!! NAME
!!
!!  hy_uhd_unsplitUpdateMultiTemp
!!
!! SYNOPSIS
!!
!!  hy_uhd_unsplitUpdateMultiTemp( integer(IN) :: blockID,
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

Subroutine hy_uhd_unsplitUpdateMultiTemp &
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
                                   hy_3Ttry_B_rad

  use Grid_interface,       ONLY : Grid_getBlkPtr, &
                                   Grid_releaseBlkPtr, &
                                   Grid_getBlkData, &
                                   Grid_getCellCoords, &
                                   Grid_getBlkIndexLimits
  
  use Logfile_interface,    ONLY : Logfile_stampMessage

  use hy_uhd_interface,     ONLY : hy_uhd_ragelike

#include "constants.h"
#include "Flash.h"
#include "Eos.h"
#include "UHD.h"

#ifdef FLASH_USM_MHD
 use Hydro_data,           ONLY : hy_hallVelocity, hy_useMagneticResistivity
#ifdef FLASH_UHD_3T
  use hy_uhd_interface,ONLY : hy_uhd_getCurrents
#endif
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
          dens_old, dens_new, &
          utot_old, utot_new, dutot_adv, &
          work_ele, work_rad, workComb, divv, Qohm  

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

  real :: lf, rf
  integer :: isize, jsize, ksize

  real    :: pele, pcombNph
  real    :: pion, pion_adv
  real    :: prad, prad_adv
  

  pcombNph = -1.0
  pion_adv = -1.0
  prad_adv = -1.0
  pele     = -1.0
  prad     = -1.0
  work_ele = -1.0
  work_rad = -1.0
  workComb = -1.0

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
              lf = (xcent(i) - dx/2)**2 / xcent(i)**2
              rf = (xcent(i) + dx/2)**2 / xcent(i)**2
           else if(hy_geometry == CYLINDRICAL) then
              lf = (xcent(i) - dx/2) / xcent(i)
              rf = (xcent(i) + dx/2) / xcent(i)
           else
              call Driver_abortFlash("[hy_uhd_unsplitUpdateMultiTemp] Unsupported Geometry")
           end if

           ! Compute change in cell internal energy over the timestep:
           dutot =  utot_new - utot_old

           ! Compute change in ion/electron/radiation internal
           ! energies due to advection:
           durad_adv = - (dt/dx)*(rf*xflux(HY_ERAD_FLUX,i+1,j,    k    ) - lf*xflux(HY_ERAD_FLUX,i,j,k)) &
                       - (dt/dy)*(   yflux(HY_ERAD_FLUX,i,  j+K2D,k    ) -    yflux(HY_ERAD_FLUX,i,j,k)) &
                       - (dt/dz)*(   zflux(HY_ERAD_FLUX,i,  j,    k+K3D) -    zflux(HY_ERAD_FLUX,i,j,k))

           duele_adv = - (dt/dx)*(rf*xflux(HY_EELE_FLUX,i+1,j,    k    ) - lf*xflux(HY_EELE_FLUX,i,j,k)) &
                       - (dt/dy)*(   yflux(HY_EELE_FLUX,i,  j+K2D,k    ) -    yflux(HY_EELE_FLUX,i,j,k)) &
                       - (dt/dz)*(   zflux(HY_EELE_FLUX,i,  j,    k+K3D) -    zflux(HY_EELE_FLUX,i,j,k))
           
           duion_adv = - (dt/dx)*(rf*xflux(HY_EION_FLUX,i+1,j,    k    ) - lf*xflux(HY_EION_FLUX,i,j,k)) &
                       - (dt/dy)*(   yflux(HY_EION_FLUX,i,  j+K2D,k    ) -    yflux(HY_EION_FLUX,i,j,k)) &
                       - (dt/dz)*(   zflux(HY_EION_FLUX,i,  j,    k+K3D) -    zflux(HY_EION_FLUX,i,j,k))

#ifdef FLASH_USM_MHD
#ifdef FLASH_UHD_3T
           if (hy_hallVelocity) then
             duele_adv = - (dt/dx)*(rf*(xflux(HY_EELE_FLUX,i+1,j,    k    ) - Jp(1,i,j,    k    )) &
                         -          lf*(xflux(HY_EELE_FLUX,i  ,j,    k    ) - Jm(1,i,j,    k    ))) &
                         - (dt/dy)*(    yflux(HY_EELE_FLUX,i,  j+K2D,k    ) - Jp(2,i,j+K2D,k    ) &
                         -              yflux(HY_EELE_FLUX,i,  j    ,k    ) - Jm(2,i,j    ,k    )) &
                         - (dt/dz)*(    zflux(HY_EELE_FLUX,i,  j,    k+K3D) - Jp(3,i,j,    k+K3D) &
                         -              zflux(HY_EELE_FLUX,i,  j,    k    ) - Jm(3,i,j,    k    ))

             durad_adv = - (dt/dx)*(rf*xflux(HY_ERAD_FLUX,i+1,j,    k    ) - lf*xflux(HY_ERAD_FLUX,i,j,k)) &
                         - (dt/dy)*(   yflux(HY_ERAD_FLUX,i,  j+K2D,k    ) -    yflux(HY_ERAD_FLUX,i,j,k)) &
                         - (dt/dz)*(   zflux(HY_ERAD_FLUX,i,  j,    k+K3D) -    zflux(HY_ERAD_FLUX,i,j,k))
           
             duion_adv = - (dt/dx)*(rf*xflux(HY_EION_FLUX,i+1,j,    k    ) - lf*xflux(HY_EION_FLUX,i,j,k)) &
                         - (dt/dy)*(   yflux(HY_EION_FLUX,i,  j+K2D,k    ) -    yflux(HY_EION_FLUX,i,j,k)) &
                         - (dt/dz)*(   zflux(HY_EION_FLUX,i,  j,    k+K3D) -    zflux(HY_EION_FLUX,i,j,k))
           endif
#endif
#endif
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
           !             - (dt/dy)*(   yflux(HY_EINT_FLUX,i,  j+K2D,k    ) -    yflux(HY_EINT_FLUX,i,j,k)) &
           !             - (dt/dz)*(   zflux(HY_EINT_FLUX,i,  j,    k+K3D) -    zflux(HY_EINT_FLUX,i,j,k))

           dutot_adv = durad_adv + duion_adv + duele_adv

           ! Compute advected ion/electron/radiation internal energies:
           uele_adv =  uele_old + duele_adv
           uion_adv =  uion_old + duion_adv
           urad_adv =  urad_old + durad_adv

           ! Check for negative advected internal energies:
           if( uele_adv <= 0.0 .or. uion_adv <= 0.0 .or. urad_adv < 0.0 ) then
              call adv_err_chk()
           end if

           !! Subtract Qohm from Eint and add it after ragelike to eele
           Qohm = 0.0              
#ifdef FLASH_USM_MHD
           if (hy_useMagneticResistivity) then
             Qohm = scrch_Ptr(HY_XN06_SCRATCHCTR_VAR,i,j,k)
           endif   
#endif

           if (hy_3TMode == HY3T_RAGELIKE) then
              call hy_uhd_ragelike(dutot-dutot_adv-Qohm, U(:,i,j,k), dens_old, &
                   duele_adv, duion_adv, durad_adv, &
                   xcent(i), ycent(j), zcent(k), &
                   uele_new, uion_new, urad_new, &
                   stepFactor=0.25*real(hy_3Ttry_G))
           else
              uele_new = uele_adv
              uion_new = uion_adv
              urad_new = urad_adv
           end if
           
           !! Update Uele with Qohm (it is actually dens*Qohm*dt at t = n)
           uele_new = uele_new + Qohm
           
           U(EELE_VAR,i,j,k) = uele_new / dens_new
           U(EION_VAR,i,j,k) = uion_new / dens_new
           U(ERAD_VAR,i,j,k) = urad_new / dens_new           

           ! The values of eion, eele, or erad may be negative or zero
           ! at this point for unknown reasons. Check for this:
           if( U(EELE_VAR, i,j,k) <= 0.0 .or. &
               U(EION_VAR, i,j,k) <= 0.0 .or. &
               U(ERAD_VAR, i,j,k) <  0.0 .or. &
               uele_adv <= 0.0 .or. &
               uion_adv <= 0.0 .or. &
               urad_adv <  0.0 ) then
               
              call Logfile_stampMessage("[hy_uhd_unsplitUpdateMultiTemp] ERROR DETECTED", .true.)
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

              call computeAddtlDiagTerms()
#if (0)
              write(errmsg, "(a,1pe20.13)") "  PELE_ADV = ", pele_adv
              call Logfile_stampMessage(trim(errmsg), .true.)
#endif

              write(errmsg, "(a,1pe20.13)") "  PION_ADV = ", pion_adv
              call Logfile_stampMessage(trim(errmsg), .true.)

#if (0)
              write(errmsg, "(a,1pe20.13)") "  PRAD_ADV = ", prad_adv
              call Logfile_stampMessage(trim(errmsg), .true.)
#endif

              write(errmsg, "(a,1pe20.13)") "  P^N+1/2  = ", pcombNph
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

              write(errmsg, "(a,1pe20.13)") "    workComb = ", workComb
              call Logfile_stampMessage(trim(errmsg), .true.)

              call Driver_abortFlash("[hy_uhd_unsplitUpdateMultiTemp] Negative 3T internal energy, CHECK LOG")
              
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

  subroutine computeAddtlDiagTerms()
    real :: avgPcomb, avgPion

    if (uele_old+duele_adv   >  0.0 .AND. &
        uion_old+duion_adv .GE. 0.0 .AND. &
        urad_old+durad_adv .GE. 0.0) then
       avgPcomb = &
                  (xflux(HY_P_FLUX,i+1,j,    k    )+xflux(HY_P_FLUX,i,j,k))
       if (NDIM>1) then
          avgPcomb = avgPcomb  &
                + (yflux(HY_P_FLUX,i,  j+K2D,k    )-yflux(HY_P_FLUX,i,j,k))
          if (NDIM>2) then
             avgPcomb = avgPcomb  &
                  + (zflux(HY_P_FLUX,i,  j,    k+K3D)-zflux(HY_P_FLUX,i,j,k))
             avgPcomb = avgPcomb  / 6.0
          else
             avgPcomb = avgPcomb  * 0.25
          end if
       else
          avgPcomb = avgPcomb  * 0.5
       end if
       avgPion = -1.0            !diag only
    else
       avgPcomb = &
                  (xflux(HY_P_FLUX,i+1,j,    k    )+xflux(HY_P_FLUX,i,j,k))
       if (NDIM>1) then
          avgPcomb = avgPcomb  &
                + (yflux(HY_P_FLUX,i,  j+K2D,k    )-yflux(HY_P_FLUX,i,j,k))
          if (NDIM>2) then
             avgPcomb = avgPcomb  &
                  + (zflux(HY_P_FLUX,i,  j,    k+K3D)-zflux(HY_P_FLUX,i,j,k))
             avgPcomb = avgPcomb  / 6.0
          else
             avgPcomb = avgPcomb  * 0.25
          end if
       else
          avgPcomb = avgPcomb  * 0.5
       end if
       avgPion = 0.0            !DEV: review later

    end if

    pcombNph = avgPcomb
    pion_adv = avgPion

    divv = &
                  (1.0/dx)* &
                 (xflux(HY_VOLU_FLUX,i+1,j,    k    )-xflux(HY_VOLU_FLUX,i,j,k)) &
                + (1.0/dy)* &
                 (yflux(HY_VOLU_FLUX,i,  j+K2D,k    )-yflux(HY_VOLU_FLUX,i,j,k)) &
                + (1.0/dz)* &
                 (zflux(HY_VOLU_FLUX,i,  j,    k+K3D)-zflux(HY_VOLU_FLUX,i,j,k))

    work_ele = &
                  (1.0/dx)*(U(PELE_VAR,i,j,k)) &
                *(xflux(HY_VOLU_FLUX,i+1,j,    k    )-xflux(HY_VOLU_FLUX,i,j,k)) &
                + (1.0/dy)*(U(PELE_VAR,i,j,k)) &
                *(yflux(HY_VOLU_FLUX,i,  j+K2D,k    )-yflux(HY_VOLU_FLUX,i,j,k)) &
                + (1.0/dz)*(U(PELE_VAR,i,j,k)) &
                *(zflux(HY_VOLU_FLUX,i,  j,    k+K3D)-zflux(HY_VOLU_FLUX,i,j,k))
    work_ele = - dt * work_ele

    work_rad = &
                  (1.0/dx)*(U(PRAD_VAR,i,j,k)) &
                *(xflux(HY_VOLU_FLUX,i+1,j,    k    )-xflux(HY_VOLU_FLUX,i,j,k)) &
                + (1.0/dy)*(U(PRAD_VAR,i,j,k)) &
                *(yflux(HY_VOLU_FLUX,i,  j+K2D,k    )-yflux(HY_VOLU_FLUX,i,j,k)) &
                + (1.0/dz)*(U(PRAD_VAR,i,j,k)) &
                *(zflux(HY_VOLU_FLUX,i,  j,    k+K3D)-zflux(HY_VOLU_FLUX,i,j,k))
    work_rad = - dt * work_rad

    workComb = &
                  (0.5/dx)*(xflux(HY_P_FLUX,i+1,j,    k    )+xflux(HY_P_FLUX,i,j,k)) &
                *(xflux(HY_VOLU_FLUX,i+1,j,    k    )-xflux(HY_VOLU_FLUX,i,j,k)) &
                + (0.5/dy)*(yflux(HY_P_FLUX,i,  j+K2D,k    )+yflux(HY_P_FLUX,i,j,k)) &
                *(yflux(HY_VOLU_FLUX,i,  j+K2D,k    )-yflux(HY_VOLU_FLUX,i,j,k)) &
                + (0.5/dz)*(zflux(HY_P_FLUX,i,  j,    k+K3D)+zflux(HY_P_FLUX,i,j,k)) &
                *(zflux(HY_VOLU_FLUX,i,  j,    k+K3D)-zflux(HY_VOLU_FLUX,i,j,k))
    workComb = - dt * workComb

    if (uele_old+duele_adv   >  0.0 .AND. &
        uion_old+duion_adv .GE. 0.0 .AND. &
        urad_old+durad_adv .GE. 0.0) then
       call hy_uhd_getPressure( &
            (uele_old+duele_adv) / dens_new, &
            (uion_old+duion_adv) / dens_new, &
            (urad_old+durad_adv) / dens_new, &
            U(:,i,j,k), &
            pele, &
            pion_adv, &
            prad)
    end if

  end subroutine computeAddtlDiagTerms


  subroutine adv_err_chk()
    implicit none

    call Logfile_stampMessage("[hy_uhd_unsplitUpdateMultiTemp] ERROR DETECTED (adv_err_chk)", .true.)
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

    call Driver_abortFlash("[hy_uhd_unsplitUpdateMultiTemp] Negative 3T internal energy, CHECK LOG")
    
  end subroutine adv_err_chk


  subroutine flux_error(flux)
    implicit none
    integer, intent(in) :: flux

    write(errmsg, "(a,2(1pe20.13))") "  Geometric Factors: ", lf, rf
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

end Subroutine hy_uhd_unsplitUpdateMultiTemp
