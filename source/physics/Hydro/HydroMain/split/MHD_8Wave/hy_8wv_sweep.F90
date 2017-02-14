!!****if* source/physics/Hydro/HydroMain/split/MHD_8Wave/hy_8wv_sweep
!!
!! NAME
!!
!!  hy_8wv_sweep
!!
!!
!! SYNOPSIS
!!
!!  hy_8wv_sweep(
!!               integer(in) :: blockCount,
!!               integer(in) :: blockList(blockCount),
!!               real(inout) :: dt,
!!               integer(in) :: sweepDir)
!!
!!
!! DESCRIPTION
!!
!!  The main 8Wave MHD routine which computes fluxes, sources
!!  and updates solution.
!!
!!
!! ARGUMENTS
!!
!!  blockCount - the number of blocks in blockList
!!  blockList  - array holding local IDs of blocks on which to advance
!!  dt         - timestep
!!  sweepDir   - direction of the sweep
!!
!!***

!!REORDER(4): U, scratchData, Rhs, totFlx



subroutine hy_8wv_sweep( blockCount, blockList, dt, sweepDir)

  use Hydro_data, ONLY : hy_cfl,  hy_xref, hy_eswitch,            &
                         hy_tref, hy_dref, hy_vref, hy_pref,      &
                         hy_eref, hy_qref, hy_bref, hy_gref,      &
                         hy_mref, hy_nref, hy_kref, hy_useGravity,&
                         hy_fluxCorrect, hy_irenorm,hy_gcMaskSize,&
                         hy_useDiffuse, hy_useMagneticResistivity,&
                         hy_useViscosity,hy_useConductivity, hy_meshMe

  use Grid_interface, ONLY : Grid_fillGuardCells,   &
                             Grid_getDeltas,        &
                             Grid_getBlkIndexLimits,&
                             Grid_getCellCoords,    &
                             Grid_getFluxData,      &
                             Grid_putFluxData,      &
                             Grid_conserveFluxes,   &
                             Grid_getBlkPtr,        &
                             Grid_releaseBlkPtr,    &
                             Grid_renormAbundance,  &
                             Grid_limitAbundance

  use Eos_interface,     ONLY : Eos_wrapped
  use Gravity_interface, ONLY : Gravity_accelOneRow

  use hy_8wv_interface, ONLY : hy_8wv_interpolate, &
                               hy_8wv_sources,     &
                               hy_8wv_fluxes,      &
                               hy_8wv_addViscousFluxes, &
                               hy_8wv_addThermalFluxes, &
                               hy_8wv_addResistiveFluxes

  use Conductivity_interface, ONLY : Conductivity
  use Viscosity_interface,    ONLY : Viscosity
  use MagneticResistivity_interface, ONLY : MagneticResistivity

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Eos.h"

  !! Argument list ---------------------------------------------
  integer, intent(IN) :: sweepDir, blockCount
  integer, intent(IN), dimension(blockCount) :: blockList
  real, intent(inout) :: dt
  !! -----------------------------------------------------------


  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC, eosRange  
  integer, dimension(MDIM)   :: dataSize
  integer, dimension(2)      :: gravPos

  integer :: sweep_eosMode = MODE_DENS_EI
  integer :: i, j, k, blockID, lnblocks, irenorm, istat, level=0
  integer :: lb, sizeX, sizeY, sizeZ, sp
  integer :: ibeg, iend, jbeg, jend, kbeg, kend
  integer :: imin, imax, jmin, jmax, kmin, kmax
  real, pointer, dimension(:,:,:,:) :: U, scratchData
  real, allocatable, dimension(:,:) :: Utemp
  real, allocatable,dimension(:)    :: xCenter, yCenter, zCenter
  real, dimension(MDIM)             :: del
  integer :: iSize, jSize, kSize, numCells
  logical :: gcell
  real    :: time, vj=0, eint, ekin, emag, dcff, smallest=TINY(1.0)
  real    :: viscKinematic

#ifdef FIXEDBLOCKSIZE
  real, dimension(NUNK_VARS,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: Rhs
  real, dimension(NFLUXES,  GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: totFlx
  real, dimension(NUNK_VARS,MAXCELLS) :: Uc, Um, Up, Src
  real, dimension(NFLUXES,  MAXCELLS) :: Flux
  real, dimension(GRID_IHI_GC)  :: dx
  real, dimension(GRID_JHI_GC)  :: dy
  real, dimension(GRID_KHI_GC)  :: dz
  real, dimension(MAXCELLS)     :: grav, speed, vint
  real, dimension(MAXCELLS)     :: viscDynamic, cond
#else
  real, allocatable,dimension(:,:,:,:) :: Rhs
  real, allocatable,dimension(:, :,:,:) :: totFlx
  real, allocatable,dimension(:,:) :: Uc, Um, Up, Src
  real, allocatable,dimension(:,:) :: Flux
  real, allocatable,dimension(:)  :: dx,dy,dz
  real, allocatable,dimension(:)  :: grav, speed, vint
  real, allocatable,dimension(:)  :: viscDynamic, cond
#endif


#ifdef FLASH_GRID_UG
  hy_fluxCorrect = .false.
#endif


  dt = dt/hy_tref
  call Grid_fillGuardCells( CENTER, ALLDIR)

  ! Main loop over leaf blocks
  do lb= 1,blockCount

     blockID = blockList(lb)
     call Grid_getDeltas(blockID, del)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     iSize=blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
     jSize=blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
     kSize=blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1
     numCells=max(iSize,jSize)
     numCells=max(numCells,kSize)
#ifndef FIXEDBLOCKSIZE
     allocate(Rhs(NUNK_VARS,iSize,jSize,kSize))
     allocate(totFlx(NFLUXES,  iSize,jSize,kSize))
     allocate(Uc(NUNK_VARS,numCells))
     allocate(Um(NUNK_VARS,numCells))
     allocate(Up(NUNK_VARS,numCells))
     allocate(Src(NUNK_VARS,numCells))
     allocate(Flux(NFLUXES,  numCells))
     allocate(dx(iSize))
     allocate(dy(jSize))
     allocate(dz(kSize))
     allocate(grav(numCells))
     allocate(speed(numCells))
     allocate(vint(numCells))
     allocate(viscDynamic(numCells))
     allocate(cond(numCells))
#endif
     if (hy_fluxCorrect) then
        ibeg=blkLimits(LOW,IAXIS)+1
        iend=blkLimits(HIGH,IAXIS)
        jbeg=blkLimits(LOW,JAXIS)+1
        jend=blkLimits(HIGH,JAXIS)
        kbeg=blkLimits(LOW,KAXIS)+1
        kend=blkLimits(HIGH,KAXIS)
     else
        ibeg=blkLimits(LOW,IAXIS)
        iend=blkLimits(HIGH,IAXIS)+1
        jbeg=blkLimits(LOW,JAXIS)
        jend=blkLimits(HIGH,JAXIS)+1
        kbeg=blkLimits(LOW,KAXIS)
        kend=blkLimits(HIGH,KAXIS)+1
     end if

     ! Set ranges of indexes
     imin=blkLimits(LOW,IAXIS)
     imax=blkLimits(HIGH,IAXIS)
     jmin=blkLimits(LOW,JAXIS)
     jmax=blkLimits(HIGH,JAXIS)
     kmin=blkLimits(LOW,KAXIS)
     kmax=blkLimits(HIGH,KAXIS)

     
     ! Setting up blkLimits to call EOS     
     eosRange = blkLimitsGC
     if(sweepDir==SWEEP_X)eosRange(HIGH,IAXIS) = blkLimits(LOW,IAXIS)-1
     if(sweepDir==SWEEP_Y)eosRange(HIGH,JAXIS) = blkLimits(LOW,JAXIS)-1
     if(sweepDir==SWEEP_Z)eosRange(HIGH,KAXIS) = blkLimits(LOW,KAXIS)-1
     call Eos_wrapped(sweep_eosMode,eosRange,blockID)
     
     eosRange = blkLimitsGC
     if(sweepDir==SWEEP_X)eosRange(LOW,IAXIS) = blkLimits(HIGH,IAXIS)+1
     if(sweepDir==SWEEP_Y)eosRange(LOW,JAXIS) = blkLimits(HIGH,JAXIS)+1
     if(sweepDir==SWEEP_Z)eosRange(LOW,KAXIS) = blkLimits(HIGH,KAXIS)+1
     call Eos_wrapped(sweep_eosMode,eosRange,blockID)
     
     
     dataSize(1) = iSize
     dataSize(2) = jSize
     dataSize(3) = kSize
     
     allocate(xCenter(iSize),stat=istat)
     allocate(yCenter(jSize),stat=istat)
     allocate(zCenter(kSize),stat=istat)
     
     gCell = .true.
     
     call Grid_getCellCoords(IAXIS, blockId, CENTER, gCell, xCenter, iSize)
     call Grid_getCellCoords(JAXIS, blockId, CENTER, gCell, yCenter, jSize)
     call Grid_getCellCoords(KAXIS, blockId, CENTER, gCell, zCenter, kSize)
     
     
     ! Get block pointer
     call Grid_getBlkPtr(blockID,U,CENTER)
     call Grid_getBlkPtr(blockID,scratchData,SCRATCH)
     
     dx(:) = del(IAXIS)
     dy(:) = del(JAXIS)
     dz(:) = del(KAXIS)

     xCenter = xCenter/hy_xref
     yCenter = yCenter/hy_xref
     zCenter = zCenter/hy_xref

     U(DENS_VAR,:,:,:) = U(DENS_VAR,:,:,:)/hy_dref
     U(ENER_VAR,:,:,:) = U(ENER_VAR,:,:,:)/hy_eref
     U(PRES_VAR,:,:,:) = U(PRES_VAR,:,:,:)/hy_pref
     U(VELX_VAR:VELZ_VAR,:,:,:) = U(VELX_VAR:VELZ_VAR,:,:,:)/hy_vref
     U(MAGX_VAR:MAGZ_VAR,:,:,:) = U(MAGX_VAR:MAGZ_VAR,:,:,:)/hy_bref


     ! Initialize arrays
     Rhs   = 0.0
     Totflx = 0.0
     grav  = 0.0
     viscDynamic  = 0.0
     cond  = 0.0
     scratchData(RESI_SCRATCH_GRID_VAR,:,:,:) = 0.0


     ! Compute right-hand side and prepare boundary fluxes
     select case(sweepDir)
     case (SWEEP_X)
        do k=kmin,kmax
           do j = jmin,jmax

              if (hy_useGravity) then
                 gravPos(1) = j
                 gravPos(2) = k
                 call Gravity_accelOneRow(gravPos,SWEEP_X,blockID,iSize,grav)
                 grav = grav/hy_gref
              endif
              
              !! store in temp array
              allocate(Utemp(NUNK_VARS,iSize),stat=istat)
              Utemp(:,:) = TRANSPOSE_IF_REORDER(U(:,:,j,k))
              call hy_8wv_interpolate(Utemp(:,:),Uc,Um,Up,grav,dt,xCenter,dx,iSize,SWEEP_X)
              deallocate(Utemp)
              call hy_8wv_sources(Uc,Um,Up,Src,grav,xCenter,dx,iSize,SWEEP_X)
              call hy_8wv_fluxes(Um,Up,Flux,speed,vint,iSize,SWEEP_X)


              if ( hy_useDiffuse ) then
                 do i = imin-1, imax+1
                    if (hy_useViscosity) then
                       call Viscosity(U(TEMP_VAR,i,j,k),U(DENS_VAR,i,j,k),&
                            U(SPECIES_BEGIN:SPECIES_END,i,j,k),viscDynamic(i),viscKinematic)
                    endif

                    if (hy_useConductivity) then
                       call Conductivity(U(TEMP_VAR,i,j,k),U(DENS_VAR,i,j,k),&
                            U(SPECIES_BEGIN:SPECIES_END,i,j,k),cond(i),dcff,2)
                    endif

                    if (hy_useMagneticResistivity) then
                       call MagneticResistivity(U(TEMP_VAR,i,j,k),U(DENS_VAR,i,j,k),&
                            U(SPECIES_BEGIN:SPECIES_END,i,j,k),scratchData(RESI_SCRATCH_GRID_VAR,i,j,k))
                    endif
                 end do

                 if (hy_useViscosity) then
                    viscDynamic    = viscDynamic/hy_nref

                    call hy_8wv_addViscousFluxes(j,k,U(VELX_VAR:VELZ_VAR,:,:,:),Flux, &
                         viscDynamic,xCenter,yCenter,zCenter,iSize,jSize,kSize,SWEEP_X)
                 endif

                 if (hy_useMagneticResistivity) then
                    scratchData(RESI_SCRATCH_GRID_VAR,:,:,:) = &
                         scratchData(RESI_SCRATCH_GRID_VAR,:,:,:)/hy_mref

                    call hy_8wv_addResistiveFluxes(j,k,U(MAGX_VAR:MAGZ_VAR,:,:,:),Flux, &
                         scratchData(RESI_SCRATCH_GRID_VAR,:,j,k),xCenter,yCenter,zCenter,&
                         iSize,jSize,kSize,SWEEP_X)
                 endif


                 if (hy_useConductivity) then
                    call hy_8wv_addThermalFluxes(j,k,U(TEMP_VAR,:,:,:),Flux, &
                         cond,xCenter,yCenter,zCenter,iSize,jSize,kSize,SWEEP_X)
                 endif


                 !! Note: Hall mhd is not supported and tested yet!
!!$                 call hy_8wv_addHallFluxes(j,k,U(DENS_VAR,:,j,k),U(MAGX_VAR:MAGZ_VAR,:,:,:), &
!!$                      Flux,vj,xCenter,yCenter,zCenter,iSize,jSize,kSize,SWEEP_X)

              end if


              Rhs(:,imin:imax,j,k) = Rhs(:,imin:imax,j,k) + TRANSPOSE_IF_REORDER(Src(:,imin:imax))

              do i = ibeg,iend
                 Rhs(DENS_VAR,i-1,j,k) = Rhs(DENS_VAR,i-1,j,k)-Flux(DENS_FLUX,i)/dx(i-1)
                 Rhs(VELX_VAR:VELZ_VAR,i-1,j,k) = &
                      Rhs(VELX_VAR:VELZ_VAR,i-1,j,k)-Flux(XMOM_FLUX:ZMOM_FLUX,i)/dx(i-1)
                 Rhs(MAGX_VAR:MAGZ_VAR,i-1,j,k) = &
                      Rhs(MAGX_VAR:MAGZ_VAR,i-1,j,k)-Flux(MAGX_FLUX:MAGZ_FLUX,i)/dx(i-1)
                 Rhs(ENER_VAR,i-1,j,k) = Rhs(ENER_VAR,i-1,j,k)-Flux(ENER_FLUX,i)/dx(i-1)

                 Rhs(DENS_VAR, i ,j,k) = Rhs(DENS_VAR, i ,j,k)+Flux(DENS_FLUX,i)/dx( i )
                 Rhs(VELX_VAR:VELZ_VAR, i ,j,k) = &
                      Rhs(VELX_VAR:VELZ_VAR, i ,j,k)+Flux(XMOM_FLUX:ZMOM_FLUX,i)/dx( i )
                 Rhs(MAGX_VAR:MAGZ_VAR, i ,j,k) = &
                      Rhs(MAGX_VAR:MAGZ_VAR, i ,j,k)+Flux(MAGX_FLUX:MAGZ_FLUX,i)/dx( i )
                 Rhs(ENER_VAR, i ,j,k) = Rhs(ENER_VAR, i ,j,k)+Flux(ENER_FLUX,i)/dx( i )
              end do

              if (hy_fluxCorrect) then
                 ! Store 1d flux array in place of 3d array
                 totFlx(1:NFLUXES,imin:imax+1,j,k)=TRANSPOSE_IF_REORDER(Flux(1:NFLUXES,imin:imax+1))
              endif

              vint = vint*dt
              do i = imin,imax
                 do sp = SPECIES_BEGIN, SPECIES_END
                    U(sp,i,j,k) = (max(0.,vint( i ))* &
                         (Up(sp,i-1)-0.5*(Up(sp,i-1)-Um(sp,i-1))* &
                         vint( i )/(dx(i-1)+vint( i )-vint(i-1)))- &
                         min(0.,vint(i+1))* &
                         (Um(sp,i+1)-0.5*(Up(sp,i+1)-Um(sp,i+1))* &
                         vint(i+1)/(dx(i+1)+vint(i+2)-vint(i+1)))+ &
                         (dx(i)+min(0.,vint(i+1))-max(0.,vint(i)))* &
                         (Uc(sp, i )-0.5*(Up(sp, i )-Um(sp, i ))* &
                         (max(0.,vint(i+1))+min(0.,vint(i)))/ &
                         (dx(i)+vint(i+1)-vint(i))))/dx(i)
                 end do

                 call hy_8wv_setTimestep(hy_meshMe,hy_cfl*hy_tref*min(dx(i)/max(speed(i),speed(i+1),vj),&
                                         0.5*dx(i)**2/max(max(viscDynamic(i),cond(i)/hy_kref)/U(DENS_VAR,i,j,k), &
                                         scratchData(RESI_SCRATCH_GRID_VAR,i,j,k),smallest)),i,j,k,blockID)
              end do

           end do
        end do
        
        if (hy_fluxCorrect) then
           call Grid_putFluxData(blockID,IAXIS,totFlx,dataSize)
        endif
        
        
#if NDIM >= 2
     case (SWEEP_Y)
        
        do k=kmin,kmax
           do i = imin,imax
              
              if (hy_useGravity) then
                 gravPos(1) = i
                 gravPos(2) = k
                 call Gravity_accelOneRow(gravPos,SWEEP_Y,blockID,jSize,grav)
                 grav = grav/hy_gref
              endif
              
              allocate(Utemp(NUNK_VARS,jSize),stat=istat)
              Utemp(:,:) = TRANSPOSE_IF_REORDER(U(:,i,:,k))
              call hy_8wv_interpolate(Utemp(:,:),Uc,Um,Up,grav,dt,yCenter,dy,jSize,SWEEP_Y)
              deallocate(Utemp)
              call hy_8wv_sources(Uc,Um,Up,Src,grav,yCenter,dy,jSize,SWEEP_Y)
              call hy_8wv_fluxes(Um,Up,Flux,speed,vint,jSize,SWEEP_Y)


              if ( hy_useDiffuse ) then
                 do j = jmin-1,jmax+1

                    if (hy_useViscosity) then
                       call Viscosity(U(TEMP_VAR,i,j,k),U(DENS_VAR,i,j,k), &
                            U(SPECIES_BEGIN:SPECIES_END,i,j,k),viscDynamic(j),viscKinematic)
                    endif

                    if (hy_useConductivity) then
                       call Conductivity(U(TEMP_VAR,i,j,k),U(DENS_VAR,i,j,k), &
                            U(SPECIES_BEGIN:SPECIES_END,i,j,k),cond(j),dcff,2)
                    endif

                    if (hy_useMagneticResistivity) then
                       call MagneticResistivity(U(TEMP_VAR,i,j,k),U(DENS_VAR,i,j,k),&
                            U(SPECIES_BEGIN:SPECIES_END,i,j,k),scratchData(RESI_SCRATCH_GRID_VAR,i,j,k))
                    endif

                 end do

                 if (hy_useViscosity) then
                    viscDynamic    = viscDynamic/hy_nref

                    call hy_8wv_addViscousFluxes(i,k,U(VELX_VAR:VELZ_VAR,:,:,:),Flux, &
                         viscDynamic,xCenter,yCenter,zCenter,iSize,jSize,kSize,SWEEP_Y)
                 endif

                 if (hy_useMagneticResistivity) then
                    scratchData(RESI_SCRATCH_GRID_VAR,:,:,:) = &
                         scratchData(RESI_SCRATCH_GRID_VAR,:,:,:)/hy_mref

                    call hy_8wv_addResistiveFluxes(i,k,U(MAGX_VAR:MAGZ_VAR,:,:,:),Flux, &
                         scratchData(RESI_SCRATCH_GRID_VAR,i,:,k),xCenter,yCenter,zCenter,&
                         iSize,jSize,kSize,SWEEP_Y)
                 endif

                 if (hy_useConductivity) then
                    call hy_8wv_addThermalFluxes(i,k,U(TEMP_VAR,:,:,:),Flux, &
                         cond,xCenter,yCenter,zCenter,iSize,jSize,kSize,SWEEP_Y)
                 endif


                 !! Hall mhd is not supported and tested yet!
!!$                 call hy_8wv_addHallFluxes(i,k,U(DENS_VAR,i,:,k),U(MAGX_VAR:MAGZ_VAR,:,:,:), &
!!$                      Flux,vj,xCenter,yCenter,zCenter,iSize,jSize,kSize,SWEEP_Y)
              end if

              Rhs(:,i,jmin:jmax,k) = Rhs(:,i,jmin:jmax,k) + TRANSPOSE_IF_REORDER(Src(:,jmin:jmax))
              
              do j=jbeg,jend
                 Rhs(DENS_VAR,i,j-1,k) = Rhs(DENS_VAR,i,j-1,k)-Flux(DENS_FLUX,j)/dy(j-1)
                 Rhs(VELX_VAR:VELZ_VAR,i,j-1,k) = &
                      Rhs(VELX_VAR:VELZ_VAR,i,j-1,k)-Flux(XMOM_FLUX:ZMOM_FLUX,j)/dy(j-1)
                 Rhs(MAGX_VAR:MAGZ_VAR,i,j-1,k) = &
                      Rhs(MAGX_VAR:MAGZ_VAR,i,j-1,k)-Flux(MAGX_FLUX:MAGZ_FLUX,j)/dy(j-1)
                 Rhs(ENER_VAR,i,j-1,k) = Rhs(ENER_VAR,i,j-1,k)-Flux(ENER_FLUX,j)/dy(j-1)

                 Rhs(DENS_VAR,i, j ,k) = Rhs(DENS_VAR,i, j ,k)+Flux(DENS_FLUX,j)/dy( j )
                 Rhs(VELX_VAR:VELZ_VAR,i, j ,k) = &
                      Rhs(VELX_VAR:VELZ_VAR,i, j ,k)+Flux(XMOM_FLUX:ZMOM_FLUX,j)/dy( j )
                 Rhs(MAGX_VAR:MAGZ_VAR,i, j ,k) = &
                      Rhs(MAGX_VAR:MAGZ_VAR,i, j ,k)+Flux(MAGX_FLUX:MAGZ_FLUX,j)/dy( j )
                 Rhs(ENER_VAR,i, j ,k) = Rhs(ENER_VAR,i, j ,k)+Flux(ENER_FLUX,j)/dy( j )
              end do

              if (hy_fluxCorrect) then
                 ! Store 1d flux array in place of 3d array
                 totFlx(1:NFLUXES,i,jmin:jmax+1,k)=TRANSPOSE_IF_REORDER(Flux(1:NFLUXES,jmin:jmax+1))
              endif

              vint = vint*dt
              do j = jmin,jmax
                 do sp = SPECIES_BEGIN, SPECIES_END
                    U(sp,i,j,k) = (max(0.,vint( j ))* &
                         (Up(sp,j-1)-0.5*(Up(sp,j-1)-Um(sp,j-1))* &
                         vint( j )/(dy(j-1)+vint( j )-vint(j-1)))- &
                         min(0.,vint(j+1))* &
                         (Um(sp,j+1)-0.5*(Up(sp,j+1)-Um(sp,j+1))* &
                         vint(j+1)/(dy(j+1)+vint(j+2)-vint(j+1)))+ &
                         (dy(j)+min(0.,vint(j+1))-max(0.,vint(j)))* &
                         (Uc(sp, j )-0.5*(Up(sp, j )-Um(sp, j ))* &
                         (max(0.,vint(j+1))+min(0.,vint(j)))/ &
                         (dy(j)+vint(j+1)-vint(j))))/dy(j)
                 end do
                 
                 ! Computing dt_mhd
                 call hy_8wv_setTimestep(hy_meshMe, hy_cfl*hy_tref*min(dy(j)/max(speed(j),speed(j+1),vj),&
                                         0.5*dy(j)**2/max(max(viscDynamic(j),cond(j)/hy_kref)/U(DENS_VAR,i,j,k), &
                                         scratchData(RESI_SCRATCH_GRID_VAR,i,j,k),smallest)),i,j,k,blockID)
              end do

           end do
        end do

        if (hy_fluxCorrect) then
           call Grid_putFluxData(blockID,JAXIS,totFlx,dataSize)
        endif

#if NDIM == 3
     case (SWEEP_Z)
        do j = jmin,jmax
           do i = imin,imax

              if (hy_useGravity) then
                 gravPos(1) = i
                 gravPos(2) = j
                 call Gravity_accelOneRow(gravPos,SWEEP_Z,blockID,kSize,grav)
                 grav = grav/hy_gref
              endif

              allocate(Utemp(NUNK_VARS,kSize),stat=istat)
              Utemp(:,:) = TRANSPOSE_IF_REORDER(U(:,i,j,:))
              call hy_8wv_interpolate(Utemp(:,:),Uc,Um,Up,grav,dt,zCenter,dz,kSize,SWEEP_Z)
              deallocate(Utemp)
              call hy_8wv_sources(Uc,Um,Up,Src,grav,zCenter,dz,kSize,SWEEP_Z)
              call hy_8wv_fluxes(Um,Up,Flux,speed,vint,kSize,SWEEP_Z)


              if ( hy_useDiffuse ) then
                 do k = kmin-1,kmax+1

                    if (hy_useViscosity) then
                       call Viscosity(U(TEMP_VAR,i,j,k),U(DENS_VAR,i,j,k), &
                            U(SPECIES_BEGIN:SPECIES_END,i,j,k),viscDynamic(k),viscKinematic)
                    endif

                    if (hy_useConductivity) then
                       call Conductivity(U(TEMP_VAR,i,j,k),U(DENS_VAR,i,j,k), &
                            U(SPECIES_BEGIN:SPECIES_END,i,j,k),cond(k),dcff,2)
                    endif

                    if (hy_useMagneticResistivity) then
                       call MagneticResistivity(U(TEMP_VAR,i,j,k),U(DENS_VAR,i,j,k),&
                            U(SPECIES_BEGIN:SPECIES_END,i,j,k),scratchData(RESI_SCRATCH_GRID_VAR,i,j,k))
                    endif

                 end do

                 if (hy_useViscosity) then
                    viscDynamic    = viscDynamic/hy_nref

                    call hy_8wv_addViscousFluxes(i,j,U(VELX_VAR:VELZ_VAR,:,:,:),Flux, &
                         viscDynamic,xCenter,yCenter,zCenter,iSize,jSize,kSize,SWEEP_Z)
                 endif

                 if (hy_useMagneticResistivity) then
                    scratchData = scratchData/hy_mref

                    call hy_8wv_addResistiveFluxes(i,j,U(MAGX_VAR:MAGZ_VAR,:,:,:),Flux, &
                         scratchData(RESI_SCRATCH_GRID_VAR,i,j,:),xCenter,yCenter,zCenter,&
                         iSize,jSize,kSize,SWEEP_Z)
                 endif

                 if (hy_useConductivity) then
                    call hy_8wv_addThermalFluxes(i,j,U(TEMP_VAR,:,:,:),Flux, &
                         cond,xCenter,yCenter,zCenter,iSize,jSize,kSize,SWEEP_Z)
                 endif


                 !! Hall mhd is not supported and tested yet!
!!$                 call hy_8wv_addHallFluxes(i,j,U(DENS_VAR,i,j,:),U(MAGX_VAR:MAGZ_VAR,:,:,:), &
!!$                      Flux,vj,xCenter,yCenter,zCenter,iSize,jSize,kSize,SWEEP_Z)

              end if

              Rhs(:,i,j,kmin:kmax) = Rhs(:,i,j,kmin:kmax) + TRANSPOSE_IF_REORDER(Src(:,kmin:kmax))
              
              do k=kbeg,kend
                 Rhs(DENS_VAR,i,j,k-1) = Rhs(DENS_VAR,i,j,k-1)-Flux(DENS_FLUX,k)/dz(k-1)
                 Rhs(VELX_VAR:VELZ_VAR,i,j,k-1) = &
                      Rhs(VELX_VAR:VELZ_VAR,i,j,k-1)-Flux(XMOM_FLUX:ZMOM_FLUX,k)/dz(k-1)
                 Rhs(MAGX_VAR:MAGZ_VAR,i,j,k-1) = &
                      Rhs(MAGX_VAR:MAGZ_VAR,i,j,k-1)-Flux(MAGX_FLUX:MAGZ_FLUX,k)/dz(k-1)
                 Rhs(ENER_VAR,i,j,k-1) = Rhs(ENER_VAR,i,j,k-1)-Flux(ENER_FLUX,k)/dz(k-1)

                 Rhs(DENS_VAR,i,j, k ) = Rhs(DENS_VAR,i,j, k )+Flux(DENS_FLUX,k)/dz( k )
                 Rhs(VELX_VAR:VELZ_VAR,i,j, k ) = &
                      Rhs(VELX_VAR:VELZ_VAR,i,j, k )+Flux(XMOM_FLUX:ZMOM_FLUX,k)/dz( k )
                 Rhs(MAGX_VAR:MAGZ_VAR,i,j, k ) = &
                      Rhs(MAGX_VAR:MAGZ_VAR,i,j, k )+Flux(MAGX_FLUX:MAGZ_FLUX,k)/dz( k )
                 Rhs(ENER_VAR,i,j, k ) = Rhs(ENER_VAR,i,j, k )+Flux(ENER_FLUX,k)/dz( k )
              end do

              if (hy_fluxCorrect) then
                 totFlx(1:NFLUXES,i,j,kmin:kmax+1)=TRANSPOSE_IF_REORDER(Flux(1:NFLUXES,kmin:kmax+1))
              endif

              vint = vint*dt
              do k = kmin,kmax
                 do sp = SPECIES_BEGIN, SPECIES_END
                    U(sp,i,j,k) = (max(0.,vint( k ))* &
                         (Up(sp,k-1)-0.5*(Up(sp,k-1)-Um(sp,k-1))* &
                         vint( k )/(dz(k-1)+vint( k )-vint(k-1)))- &
                         min(0.,vint(k+1))* &
                         (Um(sp,k+1)-0.5*(Up(sp,k+1)-Um(sp,k+1))* &
                         vint(k+1)/(dz(k+1)+vint(k+2)-vint(k+1)))+ &
                         (dz(k)+min(0.,vint(k+1))-max(0.,vint(k)))* &
                         (Uc(sp, k )-0.5*(Up(sp, k )-Um(sp, k ))* &
                         (max(0.,vint(k+1))+min(0.,vint(k)))/ &
                         (dz(k)+vint(k+1)-vint(k))))/dz(k)
                 end do
                 
                 call hy_8wv_setTimestep(hy_meshMe,hy_cfl*hy_tref*min(dz(k)/max(speed(k),speed(k+1),vj), &
                                         0.5*dz(k)**2/max(max(viscDynamic(k),cond(k)/hy_kref)/U(DENS_VAR,i,j,k), &
                                         scratchData(RESI_SCRATCH_GRID_VAR,i,j,k),smallest)),i,j,k,blockID)
              end do
           end do
        end do

        if (hy_fluxCorrect) then
           call Grid_putFluxData(blockID,KAXIS,totFlx,dataSize)
        endif
        
#endif
#endif
     end select

     ! Store internal energy
     U(EINT_VAR,:,:,:) = U(DENS_VAR,:,:,:)*(U(ENER_VAR,:,:,:)- &
          0.5*(U(VELX_VAR,:,:,:)**2+U(VELY_VAR,:,:,:)**2+U(VELZ_VAR,:,:,:)**2))
     
     ! Convert to conserved form
     U(VELX_VAR,:,:,:) = U(DENS_VAR,:,:,:)*U(VELX_VAR,:,:,:)
     U(VELY_VAR,:,:,:) = U(DENS_VAR,:,:,:)*U(VELY_VAR,:,:,:)
     U(VELZ_VAR,:,:,:) = U(DENS_VAR,:,:,:)*U(VELZ_VAR,:,:,:)
     U(ENER_VAR,:,:,:) = U(DENS_VAR,:,:,:)*U(ENER_VAR,:,:,:)+ &
          0.5*(U(MAGX_VAR,:,:,:)**2+U(MAGY_VAR,:,:,:)**2+U(MAGZ_VAR,:,:,:)**2)

     ! Update of conserved variables
     U = U+dt*Rhs

     !! Release scratch data
     call Grid_releaseBlkPtr(blockID,scratchData,SCRATCH)

     if (hy_fluxCorrect) then
        call Grid_releaseBlkPtr(blockID,U,CENTER)
     else
        ! Correct energy if necessary
        do k=kmin,kmax
           do j = jmin,jmax
              do i = imin,imax
                 ekin = 0.5*(U(VELX_VAR,i,j,k)**2+U(VELY_VAR,i,j,k)**2+U(VELZ_VAR,i,j,k)**2)&
                            /U(DENS_VAR,i,j,k)
                 emag = 0.5*(U(MAGX_VAR,i,j,k)**2+U(MAGY_VAR,i,j,k)**2+U(MAGZ_VAR,i,j,k)**2)
                 eint = U(ENER_VAR,i,j,k)-ekin-emag

                 if( eint < hy_eswitch*U(EINT_VAR,i,j,k) ) then
                    U(ENER_VAR,i,j,k) = (ekin+U(EINT_VAR,i,j,k))/U(DENS_VAR,i,j,k)
                    U(EINT_VAR,i,j,k) = U(EINT_VAR,i,j,k)/U(DENS_VAR,i,j,k)
                 else
                    U(ENER_VAR,i,j,k) = (ekin+eint)/U(DENS_VAR,i,j,k)
                    U(EINT_VAR,i,j,k) = eint/U(DENS_VAR,i,j,k)
                 end if
              end do
           end do
        end do

        ! Convert to database representation form
        U(DENS_VAR,:,:,:) = hy_dref*U(DENS_VAR,:,:,:)
        U(VELX_VAR,:,:,:) = hy_vref*U(VELX_VAR,:,:,:)/U(DENS_VAR,:,:,:)
        U(VELY_VAR,:,:,:) = hy_vref*U(VELY_VAR,:,:,:)/U(DENS_VAR,:,:,:)
        U(VELZ_VAR,:,:,:) = hy_vref*U(VELZ_VAR,:,:,:)/U(DENS_VAR,:,:,:)
        U(ENER_VAR,:,:,:) = hy_eref*U(ENER_VAR,:,:,:)
        U(MAGX_VAR:MAGZ_VAR,:,:,:) = hy_bref*U(MAGX_VAR:MAGZ_VAR,:,:,:)

        ! Renormalize or limit abundances
        if (hy_irenorm == 1) then
           call Grid_renormAbundance(blockID,blkLimits,U)
        else
           call Grid_limitAbundance(blkLimits,U)
        endif

        call Grid_releaseBlkPtr(blockID,U,CENTER)

        ! Call to Eos
        call Eos_wrapped(sweep_eosMode,blkLimits,blockID)

     endif


     ! Deallocate arrays
     deallocate(xCenter)
     deallocate(yCenter)
     deallocate(zCenter)
#ifndef FIXEDBLOCKSIZE
     deallocate(Rhs)
     deallocate(totFlx)
     deallocate(Uc)
     deallocate(Um)
     deallocate(Up)
     deallocate(Src)
     deallocate(Flux)
     deallocate(dx)
     deallocate(dy)
     deallocate(dz)
     deallocate(grav)
     deallocate(speed)
     deallocate(vint)
     deallocate(viscDynamic)
     deallocate(cond)
#endif
  end do

  if (hy_fluxCorrect) then

     ! Correct boundary fluxes
     call Grid_conserveFluxes( sweepDir, level)
     
     ! Main loop over leaf blocks
     do lb= 1,blockCount
        
        blockID = blockList(lb)
        
        call Grid_getDeltas(blockID, del)
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

        iSize=blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
        jSize=blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
        kSize=blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

        ibeg=blkLimits(LOW,IAXIS)
        iend=blkLimits(HIGH,IAXIS)
        jbeg=blkLimits(LOW,JAXIS)
        jend=blkLimits(HIGH,JAXIS)
        kbeg=blkLimits(LOW,KAXIS)
        kend=blkLimits(HIGH,KAXIS)
        
        allocate(xCenter(iSize),stat=istat)
        allocate(yCenter(jSize),stat=istat)
        allocate(zCenter(kSize),stat=istat)
#ifndef FIXEDBLOCKSIZE
     allocate(totFlx(NFLUXES,  iSize,jSize,kSize))
#endif        
        call Grid_getCellCoords(IAXIS, blockId, CENTER, gCell, xCenter, iSize)
        call Grid_getCellCoords(JAXIS, blockId, CENTER, gCell, yCenter, jSize)
        call Grid_getCellCoords(KAXIS, blockId, CENTER, gCell, zCenter, kSize)
        
        call Grid_getBlkPtr(blockID,U,CENTER)
        
        dx(:) = del(IAXIS)
        dy(:) = del(JAXIS)
        dz(:) = del(KAXIS)
        
        ! Add boundary flux contributions
        select case(sweepDir)
        case (SWEEP_X)
           call Grid_getFluxData(blockID,IAXIS,totFlx,dataSize)
           do k=kbeg,kend
              do j = jbeg,jend

                 ! Left most boundary
                 i = ibeg
                 U(DENS_VAR,i,j,k)=U(DENS_VAR,i,j,k)+(dt/dx(i))*totFlx(DENS_FLUX,i,j,k)
                 U(VELX_VAR:VELZ_VAR,i,j,k)=&
                      U(VELX_VAR:VELZ_VAR,i,j,k)+(dt/dx(i))*totFlx(XMOM_FLUX:ZMOM_FLUX,i,j,k)
                 U(MAGX_VAR:MAGZ_VAR,i,j,k)=&
                      U(MAGX_VAR:MAGZ_VAR,i,j,k)+(dt/dx(i))*totFlx(MAGX_FLUX:MAGZ_FLUX,i,j,k)
                 U(ENER_VAR,i,j,k)=U(ENER_VAR,i,j,k)+(dt/dx(i))*totFlx(ENER_FLUX,i,j,k)

                 ! Right most boundary
                 i = iend
                 U(DENS_VAR,i,j,k)=U(DENS_VAR,i,j,k)-(dt/dx(i))*totFlx(DENS_FLUX,i+1,j,k)
                 U(VELX_VAR:VELZ_VAR,i,j,k)=&
                      U(VELX_VAR:VELZ_VAR,i,j,k)-(dt/dx(i))*totFlx(XMOM_FLUX:ZMOM_FLUX,i+1,j,k)
                 U(MAGX_VAR:MAGZ_VAR,i,j,k)=&
                      U(MAGX_VAR:MAGZ_VAR,i,j,k)-(dt/dx(i))*totFlx(MAGX_FLUX:MAGZ_FLUX,i+1,j,k)
                 U(ENER_VAR,i,j,k)=U(ENER_VAR,i,j,k)-(dt/dx(i))*totFlx(ENER_FLUX,i+1,j,k)

              end do
           end do

#if NDIM >= 2
        case (SWEEP_Y)
           call Grid_getFluxData(blockID,JAXIS,totFlx,dataSize)
           do k=kbeg,kend
              do i = ibeg,iend

                 ! Left most boundary
                 j = jbeg
                 U(DENS_VAR,i,j,k)=U(DENS_VAR,i,j,k)+(dt/dy(j))*totFlx(DENS_FLUX,i,j,k)
                 U(VELX_VAR:VELZ_VAR,i,j,k)=&
                      U(VELX_VAR:VELZ_VAR,i,j,k)+(dt/dy(j))*totFlx(XMOM_FLUX:ZMOM_FLUX,i,j,k)
                 U(MAGX_VAR:MAGZ_VAR,i,j,k)=&
                      U(MAGX_VAR:MAGZ_VAR,i,j,k)+(dt/dy(j))*totFlx(MAGX_FLUX:MAGZ_FLUX,i,j,k)
                 U(ENER_VAR,i,j,k)=U(ENER_VAR,i,j,k)+(dt/dy(j))*totFlx(ENER_FLUX,i,j,k)

                 ! Right most boundary
                 j = jend
                 U(DENS_VAR,i,j,k)=U(DENS_VAR,i,j,k)-(dt/dy(j))*totFlx(DENS_FLUX,i,j+1,k)
                 U(VELX_VAR:VELZ_VAR,i,j,k)=&
                      U(VELX_VAR:VELZ_VAR,i,j,k)-(dt/dy(j))*totFlx(XMOM_FLUX:ZMOM_FLUX,i,j+1,k)
                 U(MAGX_VAR:MAGZ_VAR,i,j,k)=&
                      U(MAGX_VAR:MAGZ_VAR,i,j,k)-(dt/dy(j))*totFlx(MAGX_FLUX:MAGZ_FLUX,i,j+1,k)
                 U(ENER_VAR,i,j,k)=U(ENER_VAR,i,j,k)-(dt/dy(j))*totFlx(ENER_FLUX,i,j+1,k)

              end do
           end do
#if NDIM == 3
        case (SWEEP_Z)
           call Grid_getFluxData(blockID,KAXIS,totFlx,dataSize)
           do j = jbeg,jend
              do i = ibeg,iend

                 ! Left most boundary
                 k = kbeg
                 U(DENS_VAR,i,j,k)=U(DENS_VAR,i,j,k)+(dt/dz(k))*totFlx(DENS_FLUX,i,j,k)
                 U(VELX_VAR:VELZ_VAR,i,j,k)=&
                      U(VELX_VAR:VELZ_VAR,i,j,k)+(dt/dz(k))*totFlx(XMOM_FLUX:ZMOM_FLUX,i,j,k)
                 U(MAGX_VAR:MAGZ_VAR,i,j,k)=&
                      U(MAGX_VAR:MAGZ_VAR,i,j,k)+(dt/dz(k))*totFlx(MAGX_FLUX:MAGZ_FLUX,i,j,k)
                 U(ENER_VAR,i,j,k)=U(ENER_VAR,i,j,k)+(dt/dz(k))*totFlx(ENER_FLUX,i,j,k)

                 ! Right most boundary
                 k = kend
                 U(DENS_VAR,i,j,k)=U(DENS_VAR,i,j,k)-(dt/dz(k))*totFlx(DENS_FLUX,i,j,k+1)
                 U(VELX_VAR:VELZ_VAR,i,j,k)=&
                      U(VELX_VAR:VELZ_VAR,i,j,k)-(dt/dz(k))*totFlx(XMOM_FLUX:ZMOM_FLUX,i,j,k+1)
                 U(MAGX_VAR:MAGZ_VAR,i,j,k)=&
                      U(MAGX_VAR:MAGZ_VAR,i,j,k)-(dt/dz(k))*totFlx(MAGX_FLUX:MAGZ_FLUX,i,j,k+1)
                 U(ENER_VAR,i,j,k)=U(ENER_VAR,i,j,k)-(dt/dz(k))*totFlx(ENER_FLUX,i,j,k+1)

              end do
           end do
#endif
#endif
        end select

        ! Correct energy if necessary
        do k=kbeg,kend
           do j = jbeg,jend
              do i = ibeg,iend
                 ekin = 0.5*(U(VELX_VAR,i,j,k)**2+U(VELY_VAR,i,j,k)**2+U(VELZ_VAR,i,j,k)**2)&
                            /U(DENS_VAR,i,j,k)
                 emag = 0.5*(U(MAGX_VAR,i,j,k)**2+U(MAGY_VAR,i,j,k)**2+U(MAGZ_VAR,i,j,k)**2)
                 eint = U(ENER_VAR,i,j,k)-ekin-emag

                 if( eint < hy_eswitch*U(EINT_VAR,i,j,k) ) then
                    U(ENER_VAR,i,j,k) = (ekin+U(EINT_VAR,i,j,k))/U(DENS_VAR,i,j,k)
                    U(EINT_VAR,i,j,k) = U(EINT_VAR,i,j,k)/U(DENS_VAR,i,j,k)
                 else
                    U(ENER_VAR,i,j,k) = (ekin+eint)/U(DENS_VAR,i,j,k)
                    U(EINT_VAR,i,j,k) = eint/U(DENS_VAR,i,j,k)
                 end if
              end do
           end do
        end do

        ! Convert to database representation form
        U(DENS_VAR,:,:,:) = hy_dref*U(DENS_VAR,:,:,:)
        U(VELX_VAR,:,:,:) = hy_vref*U(VELX_VAR,:,:,:)/U(DENS_VAR,:,:,:)
        U(VELY_VAR,:,:,:) = hy_vref*U(VELY_VAR,:,:,:)/U(DENS_VAR,:,:,:)
        U(VELZ_VAR,:,:,:) = hy_vref*U(VELZ_VAR,:,:,:)/U(DENS_VAR,:,:,:)
        U(ENER_VAR,:,:,:) = hy_eref*U(ENER_VAR,:,:,:)
        U(MAGX_VAR:MAGZ_VAR,:,:,:) = hy_bref*U(MAGX_VAR:MAGZ_VAR,:,:,:)


        ! Renormalize or limit abundances
        if (hy_irenorm == 1) then
           call Grid_renormAbundance(blockID,blkLimits,U)
        else
           call Grid_limitAbundance(blkLimits,U)
        endif

        call Grid_releaseBlkPtr(blockID,U,CENTER)

        !! Call to EOS
        call Eos_wrapped(sweep_eosMode,blkLimits,blockID)

        deallocate(xCenter)
        deallocate(yCenter)
        deallocate(zCenter)
#ifndef FIXEDBLOCKSIZE
        deallocate(totFlx)
#endif       
     end do

  endif ! end of flux correction

end subroutine hy_8wv_sweep
