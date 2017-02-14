!!***if* source/physics/Hydro/explicit/split/RHD/hy_rhd_sweep
!!
!! NAME
!!
!!  hy_rhd_sweep
!!
!! SYNOPSIS
!!
!!  hy_rhd_sweep( integer (IN):: blockCount, 
!!            integer (IN):: blockList(blockCount) ,
!!            real    (IN):: timeEndAdv, 
!!            real    (IN):: dt, 
!!            real    (IN):: dtOld,            
!!            integer (IN):: sweepDir )
!!
!! DESCRIPTION
!!
!!  Performs a relativistiv hydro update of the state variables, species
!!  abundances, and mass scalars and then an update of dependant
!!  thermodynamic variables by doing a one dimensional sweep and
!!  update over a set of blocks.
!!
!! ARGUMENTS
!!
!!  blockCount -  number of blocks to advance
!!  blockList -   local block numbers of blocks to advance, sweep over these
!!  timeEndAdv -  simulation time at the end of the update
!!  dt -          timestep for advancement
!!  dtOld -       old timestep
!!  sweepDir -    direction of one-dimensional hydro sweep
!!
!!***

!!REORDER(4): U, totalFlux, Rhs

subroutine hy_rhd_sweep(blockCount,blockList,&
                     timeEndAdv, dt, dtOld, sweepDir)

  use Hydro_data, ONLY : hy_cfl,  hy_meshGeom, hy_eosMode,  &
                         hy_xref, hy_dref, hy_eref, hy_pref,&
                         hy_vref, hy_gref, hy_tref, hy_renorm, &
                         hy_fluxCorrect, hy_useGravity, hy_gamma

  use hy_rhd_interface, ONLY: hy_rhd_conserveToPrimitive, &
                              hy_rhd_primitiveToConserve, &
                              hy_rhd_states,  &
                              hy_rhd_riemann, &
                              hy_rhd_setTstep,&
                              hy_rhd_sources

  use Gravity_interface, ONLY: Gravity_accelOneRow

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

  use Eos_interface, ONLY : Eos_wrapped

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"  ! maybe not necessary
#include "RHD.h"  ! for primitive and conserved

  !! Arguments -------------------------------------------
  integer, intent(IN) ::  blockCount
  integer, intent(IN), dimension(blockCount) :: blockList
  real,    intent(IN) :: timeEndAdv, dt, dtOld
  integer, intent(IN) :: sweepDir
  !! -----------------------------------------------------


  ! Local variables
  integer :: ii,iblk,i,j,k,sp,blockID,ierr
  integer :: level=0
  integer :: iSize, jSize, kSize
  integer :: ibeg,iend,jbeg,jend,kbeg,kend
  integer, dimension(MDIM) :: dataSize
  integer, dimension(2)    :: gravPos
  integer, dimension(LOW:HIGH,MDIM) :: eosRange, blkLimits, blkLimitsGC
  real,    dimension(MDIM)    :: deltas
  real,    dimension(:,:,:,:), pointer :: U
  
  logical :: gcell = .true.


#ifdef FIXEDBLOCKSIZE
  
  real, dimension(NFLUXES,MAXCELLS) :: Flux
  real, dimension(NUNK_VARS,  MAXCELLS) :: Uc,Um,Up, Src, Utmp
  real, dimension(NFLUXES, GRID_ILO_GC:GRID_IHI_GC, &
                           GRID_JLO_GC:GRID_JHI_GC, &
                           GRID_KLO_GC:GRID_KHI_GC) :: totalFlux
  real, dimension(NUNK_VARS, GRID_ILO_GC:GRID_IHI_GC, &
                             GRID_JLO_GC:GRID_JHI_GC, &
                             GRID_KLO_GC:GRID_KHI_GC) :: Rhs
  real, dimension(GRID_ILO_GC:GRID_IHI_GC) :: xLeft, xRight, xCenter, dx, areax, dvolx
  real, dimension(GRID_JLO_GC:GRID_JHI_GC) :: yLeft, yRight, yCenter, dy, areay, dvoly
  real, dimension(GRID_KLO_GC:GRID_KHI_GC) :: zLeft, zRight, zCenter, dz, areaz, dvolz
  real, dimension(MAXCELLS) :: Dummy  ! Needed for hy_rhd_sources general interface below
  real, dimension(MAXCELLS) :: grav,speed,vint
  integer, parameter :: numCells = MAXCELLS
#else
  real, allocatable, dimension(:,:)     :: Flux,Uc,Um,Up,Src,Utmp
  real, allocatable, dimension(:,:,:,:) :: totalFlux,Rhs
  real, allocatable, dimension(:) :: xLeft,xRight,xCenter,yLeft,yRight,yCenter,zLeft,zRight,zCenter
  real, allocatable, dimension(:) :: dx, dy, dz, areax, areay, areaz, dvolx, dvoly, dvolz
  real, allocatable, dimension(:) :: Dummy
  real, allocatable, dimension(:) :: grav, speed, vint
  integer :: numCells
#endif


#ifdef FLASH_GRID_UG
  hy_fluxCorrect = .false.
#endif


  ! Fill up the guard cells to initialize the block
  ! Eos call is done below
  call Grid_fillGuardCells(CENTER,ALLDIR)

  ! Main loop over leaf blocks
  ! -----------------------------------------------------------------------
  !                   Loop on blocks begins here
  ! -----------------------------------------------------------------------
  do iblk = 1,blockCount

     blockID = blockList(iblk)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blockID, U, CENTER)
     iSize = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
     jSize = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
     kSize = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

#ifndef FIXEDBLOCKSIZE
     allocate(totalFlux(NFLUXES,iSize,jSize,kSize))
     allocate(    Rhs(NUNK_VARS,iSize,jSize,kSize))

     allocate(xLeft(iSize))
     allocate(xCenter(iSize))
     allocate(xRight(iSize))
     allocate(dx(iSize))
     allocate(areax(iSize))
     allocate(dvolx(iSize))

     allocate(yLeft(jSize))
     allocate(yCenter(jSize))
     allocate(yRight(jSize))
     allocate(dy(jSize))
     allocate(areay(jSize))
     allocate(dvoly(jSize))

     allocate(zLeft(kSize))
     allocate(zCenter(kSize))
     allocate(zRight(kSize))
     allocate(dz(kSize))
     allocate(areaz(kSize))
     allocate(dvolz(kSize))

     numCells=max(iSize,jSize)
     numCells=max(numCells,kSize)
     allocate( Utmp(NUNK_VARS, numCells))
     allocate( Uc(NUNK_VARS, numCells))
     allocate( Um(NUNK_VARS, numCells))
     allocate( Up(NUNK_VARS,numCells))
     allocate(Src(NUNK_VARS,numCells))
     allocate(Flux(NFLUXES,numCells))
     allocate(grav(numCells))
     allocate(speed(numCells))
     allocate(vint(numCells))
     allocate(Dummy(numCells))
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

     ! Running Eos on the guard cells which have just been filled
     ! Setting up blkLimits to call EOS     
     eosRange = blkLimitsGC
     if(sweepDir==SWEEP_X)eosRange(HIGH,IAXIS) = blkLimits(LOW,IAXIS)-1
     if(sweepDir==SWEEP_Y)eosRange(HIGH,JAXIS) = blkLimits(LOW,JAXIS)-1
     if(sweepDir==SWEEP_Z)eosRange(HIGH,KAXIS) = blkLimits(LOW,KAXIS)-1
     call Eos_wrapped(hy_eosMode,eosRange,blockID)
     
     eosRange = blkLimitsGC
     if(sweepDir==SWEEP_X)eosRange(LOW,IAXIS) = blkLimits(HIGH,IAXIS)+1
     if(sweepDir==SWEEP_Y)eosRange(LOW,JAXIS) = blkLimits(HIGH,JAXIS)+1
     if(sweepDir==SWEEP_Z)eosRange(LOW,KAXIS) = blkLimits(HIGH,KAXIS)+1
     call Eos_wrapped(hy_eosMode,eosRange,blockID)


     dataSize(1) = iSize
     dataSize(2) = jSize
     dataSize(3) = kSize

     ! Get coordinates
     call Grid_getCellCoords(IAXIS,blockID,LEFT_EDGE, gcell,xLeft,  iSize)
     call Grid_getCellCoords(IAXIS,blockID,CENTER,    gcell,xCenter,iSize)
     call Grid_getCellCoords(IAXIS,blockID,RIGHT_EDGE,gcell,xRight, iSize)

     call Grid_getCellCoords(JAXIS,blockID,LEFT_EDGE, gcell,yLeft,  jSize)
     call Grid_getCellCoords(JAXIS,blockID,CENTER,    gcell,yCenter,jSize)
     call Grid_getCellCoords(JAXIS,blockID,RIGHT_EDGE,gcell,yRight, jSize)

     call Grid_getCellCoords(KAXIS,blockID,LEFT_EDGE, gcell,zLeft,  kSize)
     call Grid_getCellCoords(KAXIS,blockID,CENTER,    gcell,zCenter,kSize)
     call Grid_getCellCoords(KAXIS,blockID,RIGHT_EDGE,gcell,zRight, kSize)


     ! Get block data
     call Grid_getBlkPtr(blockID,U,CENTER)

     ! Get deltas
     call Grid_getDeltas(blockID,deltas)
     dx(:) = deltas(IAXIS)
     dy(:) = deltas(JAXIS)
     dz(:) = deltas(KAXIS)

     if (hy_meshGeom == CARTESIAN) then 
        dvolx = dx
        dvoly = dy
        dvolz = dz
        areax = 1.0
        areay = 1.0
        areaz = 1.0
     else if (hy_meshGeom == CYLINDRICAL) then
        dvolx = xCenter*dx ! r*dr
        dvoly = dy
        dvolz = 0.0   ! CYLINDRICAL only works in 2-D
        areax = xLeft
        areay = 1.0
        areaz = 0.0   ! CYLINDRICAL only works in 2-D 
     else if (hy_meshGeom == SPHERICAL) then
        dvolx = ((xCenter + 0.5*dx)**3.0 - (xCenter - 0.5*dx)**3.0)/3.0 
        dvoly = 0.0   ! SPHERICAL only works in 1-D
        dvolz = 0.0   ! SPHERICAL only works in 1-D
        areax = xLeft*xLeft
        areay = 0.0   ! SPHERICAL only works in 1-D
        areaz = 0.0   ! SPHERICAL only works in 1-D
     end if

     ! Scaling
     xLeft  = xLeft /hy_xref
     xRight = xRight/hy_xref
     yLeft  = yLeft /hy_xref
     yRight = yRight/hy_xref
     zLeft  = zLeft /hy_xref
     zRight = zRight/hy_xref

     ! Scaling by relativistic values
     U(DENS_VAR,:,:,:) = U(DENS_VAR,:,:,:)/hy_dref
     U(ENER_VAR,:,:,:) = U(ENER_VAR,:,:,:)/hy_eref
     U(PRES_VAR,:,:,:) = U(PRES_VAR,:,:,:)/hy_pref
     U(VELX_VAR:VELZ_VAR,:,:,:) = U(VELX_VAR:VELZ_VAR,:,:,:)/hy_vref

     ! Initialize to zero
     Rhs   = 0.0
     grav  = 0.0
     Src   = 0.0
     
     ! Compute right-hand side and prepare boundary fluxes
     select case(sweepDir)
     case (SWEEP_X)
        do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)

              if (hy_useGravity) then
                 gravPos(1) = j
                 gravPos(2) = k
                 call Gravity_accelOneRow(gravPos,SWEEP_X,blockID,iSize,grav)
                 grav = grav/hy_gref
              endif

              ! One dimensional slice
              Uc(:,:) = TRANSPOSE_IF_REORDER(U(:,:,j,k))

              call hy_rhd_states (Uc, Um, Up, grav, dt, xCenter, dx, dvolx, iSize, SWEEP_X)
              call hy_rhd_riemann(Uc, Um, Up, Flux, xCenter, dx, dt, speed, vint, iSize, SWEEP_X)
              call hy_rhd_sources(Uc, Dummy, Dummy, Src, xCenter, dx, &
                               blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                               iSize, SWEEP_X, HY_CONSERVE)

              Rhs(:,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),j,k) = &
                   Rhs(:,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),j,k)+&
                   TRANSPOSE_IF_REORDER(Src(:,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)))

              if (hy_meshGeom .NE. CARTESIAN) then
                 ! -- multiply fluxes with area  ----------------- !
                 do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)+1
                    Flux(:,i) = Flux(:,i)*areax(i)
                 end do
              endif

              do i = ibeg,iend
                 Rhs(DENS_VAR,i-1,j,k) = Rhs(DENS_VAR,i-1,j,k)-Flux(DENS_FLUX,i)/dvolx(i-1)
                 Rhs(DENS_VAR,i  ,j,k) = Rhs(DENS_VAR,i  ,j,k)+Flux(DENS_FLUX,i)/dvolx(i  )

                 Rhs(VELX_VAR,i-1,j,k) = Rhs(VELX_VAR,i-1,j,k)-Flux(XMOM_FLUX,i)/dvolx(i-1)
                 Rhs(VELX_VAR,i  ,j,k) = Rhs(VELX_VAR,i  ,j,k)+Flux(XMOM_FLUX,i)/dvolx(i  )

                 Rhs(VELY_VAR,i-1,j,k) = Rhs(VELY_VAR,i-1,j,k)-Flux(YMOM_FLUX,i)/dvolx(i-1)
                 Rhs(VELY_VAR,i  ,j,k) = Rhs(VELY_VAR,i  ,j,k)+Flux(YMOM_FLUX,i)/dvolx(i  )

                 Rhs(VELZ_VAR,i-1,j,k) = Rhs(VELZ_VAR,i-1,j,k)-Flux(ZMOM_FLUX,i)/dvolx(i-1)
                 Rhs(VELZ_VAR,i  ,j,k) = Rhs(VELZ_VAR,i  ,j,k)+Flux(ZMOM_FLUX,i)/dvolx(i  )

                 Rhs(ENER_VAR,i-1,j,k) = Rhs(ENER_VAR,i-1,j,k)-Flux(ENER_FLUX,i)/dvolx(i-1)
                 Rhs(ENER_VAR,i  ,j,k) = Rhs(ENER_VAR,i  ,j,k)+Flux(ENER_FLUX,i)/dvolx(i  )
              end do

              if (hy_fluxCorrect) then
                 !! Store the flux for later correction
                 totalFlux(1:NFLUXES,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+1,j,k)=&
                      TRANSPOSE_IF_REORDER(Flux(1:NFLUXES,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+1))
              endif

              ! Update the interface velocity with the time step
              vint = vint*dt
              do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
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
                 call hy_rhd_setTstep(hy_cfl*hy_tref*dx(i)/max(speed(i),speed(i+1)),i,j,k,blockID)

              end do
           end do
        end do

        if (hy_fluxCorrect) then
           call Grid_putFluxData(blockID,IAXIS,totalFlux,dataSize)
        endif

#if NDIM >= 2
     case (SWEEP_Y)
        do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
           do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

              if (hy_useGravity) then
                 gravPos(1) = i
                 gravPos(2) = k
                 call Gravity_accelOneRow(gravPos,SWEEP_Y,blockID,jSize,grav)
                 grav = grav/hy_gref
              endif

              ! One dimensional slice
              Uc(:,:) = TRANSPOSE_IF_REORDER(U(:,i,:,k))
              call hy_rhd_states (Uc, Um, Up, grav, dt, yCenter, dy, dvoly, jSize, SWEEP_Y)
              call hy_rhd_riemann(Uc, Um, Up, Flux, yCenter, dy, dt, speed, vint, jSize, SWEEP_Y)
              call hy_rhd_sources(Uc, Dummy, Dummy, Src, yCenter, dy, &
                               blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                               jSize, SWEEP_Y, HY_CONSERVE)

              Rhs(:,i,blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),k) = &
                   Rhs(:,i,blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),k)+&
                   TRANSPOSE_IF_REORDER(Src(:,blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)))

              if (hy_meshGeom .NE. CARTESIAN) then
                 ! -- multiply fluxes with area  ----------------- !
                 do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)+1
                    Flux(:,j) = Flux(:,j)*areay(j)
                 end do
              endif

              do j = jbeg,jend
                 Rhs(DENS_VAR,i,j-1,k) = Rhs(DENS_VAR,i,j-1,k)-Flux(DENS_FLUX,j)/dvoly(j-1)
                 Rhs(DENS_VAR,i, j ,k) = Rhs(DENS_VAR,i, j ,k)+Flux(DENS_FLUX,j)/dvoly( j )

                 Rhs(VELX_VAR,i,j-1,k) = Rhs(VELX_VAR,i,j-1,k)-Flux(XMOM_FLUX,j)/dvoly(j-1)
                 Rhs(VELX_VAR,i, j ,k) = Rhs(VELX_VAR,i, j ,k)+Flux(XMOM_FLUX,j)/dvoly( j )

                 Rhs(VELY_VAR,i,j-1,k) = Rhs(VELY_VAR,i,j-1,k)-Flux(YMOM_FLUX,j)/dvoly(j-1)
                 Rhs(VELY_VAR,i, j ,k) = Rhs(VELY_VAR,i, j ,k)+Flux(YMOM_FLUX,j)/dvoly( j )

                 Rhs(VELZ_VAR,i,j-1,k) = Rhs(VELZ_VAR,i,j-1,k)-Flux(ZMOM_FLUX,j)/dvoly(j-1)
                 Rhs(VELZ_VAR,i, j ,k) = Rhs(VELZ_VAR,i, j ,k)+Flux(ZMOM_FLUX,j)/dvoly( j )

                 Rhs(ENER_VAR,i,j-1,k) = Rhs(ENER_VAR,i,j-1,k)-Flux(ENER_FLUX,j)/dvoly(j-1)
                 Rhs(ENER_VAR,i, j ,k) = Rhs(ENER_VAR,i, j ,k)+Flux(ENER_FLUX,j)/dvoly( j )
              end do

              if (hy_fluxCorrect) then
                 !! Store the flux for later correction
                 totalFlux(1:NFLUXES,i,blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+1,k)=&
                      TRANSPOSE_IF_REORDER(Flux(1:NFLUXES,blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+1))
              endif

              ! Update the interface velocity with the time step
              vint = vint*dt
              do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
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

                 call hy_rhd_setTstep(hy_cfl*hy_tref*dy(j)/max(speed(j),speed(j+1)),i,j,k,blockID)

              end do
           end do
        end do

        if (hy_fluxCorrect) then
           call Grid_putFluxData(blockID,JAXIS,totalFlux,dataSize)
        endif

#if NDIM == 3
     case (SWEEP_Z)
        do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

              if (hy_useGravity) then
                 gravPos(1) = i
                 gravPos(2) = j
                 call Gravity_accelOneRow(gravPos,SWEEP_Z,blockID,kSize,grav)
                 grav = grav/hy_gref
              endif

              ! One dimensional slice
              Uc(:,:) = TRANSPOSE_IF_REORDER(U(:,i,j,:))
              call hy_rhd_states (Uc, Um, Up, grav, dt, zCenter, dz, dvolz, kSize, SWEEP_Z)
              call hy_rhd_riemann(Uc, Um, Up, Flux, zCenter, dz, dt, speed, vint, kSize, SWEEP_Z)
              call hy_rhd_sources(Uc, Dummy, Dummy, Src, zCenter, dz, &
                               blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),&
                               kSize, SWEEP_Z, HY_CONSERVE)

              Rhs(:,i,j,blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = &
                   Rhs(:,i,j,blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))+&
                   TRANSPOSE_IF_REORDER(Src(:,blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)))

              if (hy_meshGeom .NE. CARTESIAN) then
                 ! -- multiply fluxes with area  ----------------- !
                 do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)+1
                    Flux(:,k) = Flux(:,k)*areaz(k)
                 end do
              endif

              do k = kbeg,kend
                 Rhs(DENS_VAR,i,j,k-1) = Rhs(DENS_VAR,i,j,k-1)-Flux(DENS_FLUX,k)/dvolz(k-1)
                 Rhs(DENS_VAR,i,j, k ) = Rhs(DENS_VAR,i,j, k )+Flux(DENS_FLUX,k)/dvolz( k )

                 Rhs(VELX_VAR,i,j,k-1) = Rhs(VELX_VAR,i,j,k-1)-Flux(XMOM_FLUX,k)/dvolz(k-1)
                 Rhs(VELX_VAR,i,j, k ) = Rhs(VELX_VAR,i,j, k )+Flux(XMOM_FLUX,k)/dvolz( k )

                 Rhs(VELY_VAR,i,j,k-1) = Rhs(VELY_VAR,i,j,k-1)-Flux(YMOM_FLUX,k)/dvolz(k-1)
                 Rhs(VELY_VAR,i,j, k ) = Rhs(VELY_VAR,i,j, k )+Flux(YMOM_FLUX,k)/dvolz( k )

                 Rhs(VELZ_VAR,i,j,k-1) = Rhs(VELZ_VAR,i,j,k-1)-Flux(ZMOM_FLUX,k)/dvolz(k-1)
                 Rhs(VELZ_VAR,i,j, k ) = Rhs(VELZ_VAR,i,j, k )+Flux(ZMOM_FLUX,k)/dvolz( k )

                 Rhs(ENER_VAR,i,j,k-1) = Rhs(ENER_VAR,i,j,k-1)-Flux(ENER_FLUX,k)/dvolz(k-1)
                 Rhs(ENER_VAR,i,j, k ) = Rhs(ENER_VAR,i,j, k )+Flux(ENER_FLUX,k)/dvolz( k )
              end do

              if (hy_fluxCorrect) then
                 !! Store the flux for later correction
                 totalFlux(1:NFLUXES,i,j,blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)+1)=&
                      TRANSPOSE_IF_REORDER(Flux(1:NFLUXES,blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)+1)))
              endif

              ! Update the interface velocity with the time step
              vint = vint*dt
              do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
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

                 call hy_rhd_setTstep(hy_cfl*hy_tref*dz(k)/max(speed(k),speed(k+1)),i,j,k,blockID)

              end do
           end do
        end do

        if (hy_fluxCorrect) then
           call Grid_putFluxData(blockID,KAXIS,totalFlux,dataSize)
        endif

#endif
#endif
     end select


     ! Convert to conserved form before proceeding temporal update
     select case(sweepDir)
     case (SWEEP_X)
        do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
              Utmp(:,:) = TRANSPOSE_IF_REORDER(U(:,:,j,k))
              call hy_rhd_primitiveToConserve (Utmp,blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),iSize)
              U(:,:,j,k) = TRANSPOSE_IF_REORDER(Utmp(:,:))
           end do
        end do
     case (SWEEP_Y)
        do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
           do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              Utmp(:,:) = TRANSPOSE_IF_REORDER(U(:,i,:,k))
              call hy_rhd_primitiveToConserve (Utmp,blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS), jSize)
              U(:,i,:,k) = TRANSPOSE_IF_REORDER(Utmp(:,:))
           end do
        end do
     case (SWEEP_Z)
        do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              Utmp(:,:) = TRANSPOSE_IF_REORDER(U(:,i,j,:))
              call hy_rhd_primitiveToConserve (Utmp,blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS), kSize)
              U(:,i,j,:) = TRANSPOSE_IF_REORDER(Utmp(:,:))
           end do
        end do
     end select


     !-----------------------------------------------------------------------
     ! Update conserved variables to the next time step

     U = U + dt * Rhs

     !-----------------------------------------------------------------------

     ! Release pointer here and will do flux correction later
     if (hy_fluxCorrect) then
        call Grid_releaseBlkPtr(blockID,U,CENTER)
     else

        ! if no flux correction is needed then we finalize calculation here
        ! Convert conserved form to primitive form
        select case(sweepDir)
        case (SWEEP_X)
           do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
              do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)  
                 Utmp(:,:) = TRANSPOSE_IF_REORDER(U(:,:,j,k))
                 call hy_rhd_conserveToPrimitive (Utmp,blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS), iSize)
                 U(:,:,j,k) = TRANSPOSE_IF_REORDER(Utmp(:,:))
              end do
           end do
        case (SWEEP_Y)
           do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
              do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
                 Utmp(:,:) = TRANSPOSE_IF_REORDER(U(:,i,:,k))
                 call hy_rhd_conserveToPrimitive (Utmp,blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS), jSize)
                 U(:,i,:,k) = TRANSPOSE_IF_REORDER(Utmp(:,:))
              end do
           end do
        case (SWEEP_Z)
           do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
                 Utmp(:,:) = TRANSPOSE_IF_REORDER(U(:,i,j,:))
                 call hy_rhd_conserveToPrimitive (Utmp,blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS), kSize) 
                 U(:,i,j,:) = TRANSPOSE_IF_REORDER(Utmp(:,:))
              end do
           end do
        end select

        ! Note: internal energy is taken care of in RHD eos routine
        U(DENS_VAR,:,:,:) = hy_dref*U(DENS_VAR,:,:,:)
        U(VELX_VAR:VELZ_VAR,:,:,:) = hy_vref*U(VELX_VAR:VELZ_VAR,:,:,:)
        U(PRES_VAR,:,:,:) = hy_pref*U(PRES_VAR,:,:,:)
        U(ENER_VAR,:,:,:) = hy_eref*U(ENER_VAR,:,:,:)/U(DENS_VAR,:,:,:)

        ! Release pointer and finalize calculation
        call Grid_releaseBlkPtr(blockID,U,CENTER)

        ! Call to Eos before closing
        call Eos_wrapped(hy_eosMode,blkLimits,blockID)
     endif

#ifndef FIXEDBLOCKSIZE
     deallocate(totalFlux)
     deallocate(Rhs)

     deallocate(xLeft)
     deallocate(xCenter)
     deallocate(xRight)
     deallocate(dx)
     deallocate(areax)
     deallocate(dvolx)

     deallocate(yLeft)
     deallocate(yCenter)
     deallocate(yRight)
     deallocate(dy)
     deallocate(areay)
     deallocate(dvoly)

     deallocate(zLeft)
     deallocate(zCenter)
     deallocate(zRight)
     deallocate(dz)
     deallocate(areaz)
     deallocate(dvolz)

     deallocate(Utmp)
     deallocate(Uc)
     deallocate(Um)
     deallocate(Up)
     deallocate(Src)
     deallocate(Flux)
     deallocate(grav)
     deallocate(speed)
     deallocate(vint)
     deallocate(Dummy)
#endif
  end do
  ! ---------------------------------------------------------------------
  !                   Loop on blocks ends here  
  ! ---------------------------------------------------------------------


  ! ---------------------------------------------------------------------
  !     Now begin the second part to correct the fluxes at the edges
  ! ---------------------------------------------------------------------
  if (hy_fluxCorrect) then
     call Grid_conserveFluxes(sweepDir,level)

     ! Main loop over leaf blocks
     do iblk = 1,blockCount

        blockID=blockList(iblk)

        call Grid_getDeltas(blockID,deltas)
        dx = deltas(IAXIS)
        dy = deltas(JAXIS)
        dz = deltas(KAXIS)

        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
        iSize = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
        jSize = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
        kSize = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

        ibeg=blkLimits(LOW,IAXIS)
        iend=blkLimits(HIGH,IAXIS)
        jbeg=blkLimits(LOW,JAXIS)
        jend=blkLimits(HIGH,JAXIS)
        kbeg=blkLimits(LOW,KAXIS)
        kend=blkLimits(HIGH,KAXIS)


#ifndef FIXEDBLOCKSIZE
        allocate(totalFlux(NFLUXES,iSize,jSize,kSize))

        allocate(xCenter(iSize))
        allocate(dx(iSize))
        allocate(dvolx(iSize))

        allocate(yCenter(jSize))
        allocate(dy(jSize))
        allocate(dvoly(jSize))

        allocate(zCenter(kSize))
        allocate(dz(kSize))
        allocate(dvolz(kSize))
#endif

        ! Get coordinates
        call Grid_getCellCoords(IAXIS,blockID,CENTER,gcell,xCenter,iSize)
        call Grid_getCellCoords(JAXIS,blockID,CENTER,gcell,yCenter,jSize)
        call Grid_getCellCoords(KAXIS,blockID,CENTER,gcell,zCenter,kSize)

        ! The original version only considered changes in IAXIS
        ! because that is the only one that differs from dx
        ! but here we repeat it all to avoid confusion
        if (hy_meshGeom == CARTESIAN) then 
           dvolx = dx
           dvoly = dy
           dvolz = dz
        else if (hy_meshGeom == CYLINDRICAL) then
           dvolx = xCenter*dx
           dvoly = dy
           dvolz = 0.0   ! CYLINDRICAL only works in 2-D
        else if (hy_meshGeom == SPHERICAL) then
           dvolx = ((xCenter + 0.5*dx)**3.0 - (xCenter - 0.5*dx)**3.0)/3.0
           dvoly = 0.0   ! SPHERICAL only works in 1-D
           dvolz = 0.0   ! SPHERICAL only works in 1-D
        end if

        ! Add boundary flux contributions
        dataSize(IAXIS)=iSize
        dataSize(JAXIS)=jSize
        dataSize(KAXIS)=kSize

        call Grid_getFluxData(blockID,sweepDir,totalFlux,dataSize)
        call Grid_getBlkPtr(blockID,U,CENTER)

        select case(sweepDir)
        case (SWEEP_X)
           ! Do left and right boundaries
           do k=kbeg,kend
              do j = jbeg,jend
                 U(DENS_VAR,ibeg,j,k) = U(DENS_VAR,ibeg,j,k) + dt/dvolx(ibeg)*totalFlux(DENS_FLUX,ibeg  ,j,k)
                 U(DENS_VAR,iend,j,k) = U(DENS_VAR,iend,j,k) - dt/dvolx(iend)*totalFlux(DENS_FLUX,iend+1,j,k)

                 U(VELX_VAR,ibeg,j,k) = U(VELX_VAR,ibeg,j,k) + dt/dvolx(ibeg)*totalFlux(XMOM_FLUX,ibeg,  j,k)
                 U(VELX_VAR,iend,j,k) = U(VELX_VAR,iend,j,k) - dt/dvolx(iend)*totalFlux(XMOM_FLUX,iend+1,j,k)

                 U(VELY_VAR,ibeg,j,k) = U(VELY_VAR,ibeg,j,k) + dt/dvolx(ibeg)*totalFlux(YMOM_FLUX,ibeg,  j,k)
                 U(VELY_VAR,iend,j,k) = U(VELY_VAR,iend,j,k) - dt/dvolx(iend)*totalFlux(YMOM_FLUX,iend+1,j,k)

                 U(VELZ_VAR,ibeg,j,k) = U(VELZ_VAR,ibeg,j,k) + dt/dvolx(ibeg)*totalFlux(ZMOM_FLUX,ibeg,  j,k)
                 U(VELZ_VAR,iend,j,k) = U(VELZ_VAR,iend,j,k) - dt/dvolx(iend)*totalFlux(ZMOM_FLUX,iend+1,j,k)

                 U(ENER_VAR,ibeg,j,k) = U(ENER_VAR,ibeg,j,k) + dt/dvolx(ibeg)*totalFlux(ENER_FLUX,ibeg,  j,k)
                 U(ENER_VAR,iend,j,k) = U(ENER_VAR,iend,j,k) - dt/dvolx(iend)*totalFlux(ENER_FLUX,iend+1,j,k)
              end do
           end do

#if NDIM >= 2
        case (SWEEP_Y)
           do k=kbeg,kend
              do i = ibeg,iend
                 U(DENS_VAR,i,jbeg,k) = U(DENS_VAR,i,jbeg,k) + dt/dvoly(jbeg)*totalFlux(DENS_FLUX,i,jbeg  ,k)
                 U(DENS_VAR,i,jend,k) = U(DENS_VAR,i,jend,k) - dt/dvoly(jend)*totalFlux(DENS_FLUX,i,jend+1,k)

                 U(VELX_VAR,i,jbeg,k) = U(VELX_VAR,i,jbeg,k) + dt/dvoly(jbeg)*totalFlux(XMOM_FLUX,i,jbeg  ,k)
                 U(VELX_VAR,i,jend,k) = U(VELX_VAR,i,jend,k) - dt/dvoly(jend)*totalFlux(XMOM_FLUX,i,jend+1,k)

                 U(VELY_VAR,i,jbeg,k) = U(VELY_VAR,i,jbeg,k) + dt/dvoly(jbeg)*totalFlux(YMOM_FLUX,i,jbeg  ,k)
                 U(VELY_VAR,i,jend,k) = U(VELY_VAR,i,jend,k) - dt/dvoly(jend)*totalFlux(YMOM_FLUX,i,jend+1,k)

                 U(VELZ_VAR,i,jbeg,k) = U(VELZ_VAR,i,jbeg,k) + dt/dvoly(jbeg)*totalFlux(ZMOM_FLUX,i,jbeg  ,k)
                 U(VELZ_VAR,i,jend,k) = U(VELZ_VAR,i,jend,k) - dt/dvoly(jend)*totalFlux(ZMOM_FLUX,i,jend+1,k)

                 U(ENER_VAR,i,jbeg,k) = U(ENER_VAR,i,jbeg,k) + dt/dvoly(jbeg)*totalFlux(ENER_FLUX,i,jbeg  ,k)
                 U(ENER_VAR,i,jend,k) = U(ENER_VAR,i,jend,k) - dt/dvoly(jend)*totalFlux(ENER_FLUX,i,jend+1,k)
              end do
           end do

#if NDIM == 3
        case (SWEEP_Z)
           do j = jbeg,jend
              do i = ibeg,iend
                 U(DENS_VAR,i,j,kbeg) = U(DENS_VAR,i,j,kbeg) + dt/dvolz(kbeg)*totalFlux(DENS_FLUX,i,j,kbeg  )
                 U(DENS_VAR,i,j,kend) = U(DENS_VAR,i,j,kend) - dt/dvolz(kend)*totalFlux(DENS_FLUX,i,j,kend+1)

                 U(VELX_VAR,i,j,kbeg) = U(VELX_VAR,i,j,kbeg) + dt/dvolz(kbeg)*totalFlux(XMOM_FLUX,i,j,kbeg  )
                 U(VELX_VAR,i,j,kend) = U(VELX_VAR,i,j,kend) - dt/dvolz(kend)*totalFlux(XMOM_FLUX,i,j,kend+1)

                 U(VELY_VAR,i,j,kbeg) = U(VELY_VAR,i,j,kbeg) + dt/dvolz(kbeg)*totalFlux(YMOM_FLUX,i,j,kbeg  )
                 U(VELY_VAR,i,j,kend) = U(VELY_VAR,i,j,kend) - dt/dvolz(kend)*totalFlux(YMOM_FLUX,i,j,kend+1)

                 U(VELZ_VAR,i,j,kbeg) = U(VELZ_VAR,i,j,kbeg) + dt/dvolz(kbeg)*totalFlux(ZMOM_FLUX,i,j,kbeg  )
                 U(VELZ_VAR,i,j,kend) = U(VELZ_VAR,i,j,kend) - dt/dvolz(kend)*totalFlux(ZMOM_FLUX,i,j,kend+1)

                 U(ENER_VAR,i,j,kbeg) = U(ENER_VAR,i,j,kbeg) + dt/dvolz(kbeg)*totalFlux(ENER_FLUX,i,j,kbeg  )
                 U(ENER_VAR,i,j,kend) = U(ENER_VAR,i,j,kend) - dt/dvolz(kend)*totalFlux(ENER_FLUX,i,j,kend+1)
              end do
           end do
#endif
#endif
        end select

        select case(sweepDir)
        case (SWEEP_X)
           do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
              do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
                 Utmp(:,:) = TRANSPOSE_IF_REORDER(U(:,:,j,k))
                 call hy_rhd_conserveToPrimitive (Utmp,blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS), iSize)
                 U(:,:,j,k) = TRANSPOSE_IF_REORDER(Utmp(:,:))
              end do
           end do
        case (SWEEP_Y)
           do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
              do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
                 Utmp(:,:) = TRANSPOSE_IF_REORDER(U(:,i,:,k))
                 call hy_rhd_conserveToPrimitive (Utmp,blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS), jSize)
                 U(:,i,:,k) = TRANSPOSE_IF_REORDER(Utmp(:,:))
              end do
           end do
        case (SWEEP_Z)
           do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
                 Utmp(:,:) = TRANSPOSE_IF_REORDER(U(:,i,j,:))
                 call hy_rhd_conserveToPrimitive (Utmp,blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS), kSize)
                 U(:,i,j,:) = TRANSPOSE_IF_REORDER(Utmp(:,:))
              end do
           end do
        end select

        ! Note: internal energy is taken care of in RHD eos routine
        U(DENS_VAR,:,:,:) = hy_dref*U(DENS_VAR,:,:,:)
        U(VELX_VAR:VELZ_VAR,:,:,:) = hy_vref*U(VELX_VAR:VELZ_VAR,:,:,:)
        U(PRES_VAR,:,:,:) = hy_pref*U(PRES_VAR,:,:,:)
        U(ENER_VAR,:,:,:) = hy_eref*U(ENER_VAR,:,:,:)/U(DENS_VAR,:,:,:)

    
#ifndef FIXEDBLOCKSIZE
        deallocate(totalFlux)

        deallocate(xCenter)
        deallocate(dx)
        deallocate(dvolx)

        deallocate(yCenter)
        deallocate(dy)
        deallocate(dvoly)

        deallocate(zCenter)
        deallocate(dz)
        deallocate(dvolz)
#endif

        ! Renormalize or limit abundances
        if (hy_renorm == 1) then
           call Grid_renormAbundance(blockID,blkLimits,U)
        else
           call Grid_limitAbundance(blkLimits,U)
        endif

        ! Release Block pointer
        call Grid_releaseBlkPtr(blockID,U,CENTER)
 
        ! Call to Eos
        call Eos_wrapped(hy_eosMode,blkLimits,blockID)

     end do ! end of flux correction over blockID
  end if

end subroutine hy_rhd_sweep
