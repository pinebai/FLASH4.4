!!****if* source/physics/sourceTerms/Stir/StirMain/FromFile/Stir
!!
!! NAME
!!  Stir
!!
!! SYNOPSIS
!!  Stir(integer(IN) :: blockCount,
!!       integer(IN) :: blockList(blockCount),
!!       real(IN)    :: dt)
!!
!! DESCRIPTION
!!   Apply the stirring operator on the list of blocks provided as input
!!
!! ARGUMENTS
!!   blockCount   : The number of blocks in the list
!!   blockList(:) : The list of blocks on which to apply the stirring operator
!!   dt           : the current timestep
!!
!! AUTHOR
!!   Christoph Federrath, 2008-2013
!!
!!    (Aug 2013: added a write-out of the correction for the center-of-mass motion)
!!    (2012: added a pre-proc statement to treat the acceleration field as a force for testing; not the default)
!!    (2011/2012: added a write-out of the energy injection rate)
!!
!!***

subroutine Stir(blockCount,blockList,dt)

  use Stir_data,        ONLY : st_useStir, st_infilename, st_lastTimeUpdatedAccel, &
                               st_dtUpdateAccel, &
                               st_globalMe, st_meshComm
  use Driver_interface, ONLY : Driver_getSimTime
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_interface,   ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, &
                               Grid_getDeltas, Grid_releaseBlkPtr, Grid_getCellCoords

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer, intent(IN)                        :: blockCount
  integer, dimension(blockCount), intent(IN) :: blockList
  real, intent(IN)                           :: dt

  logical, save                :: firstCall = .true.
  integer                      :: blockID, i, j, k
  integer, dimension(2,MDIM)   :: blkLimits, blkLimitsGC
  logical, parameter           :: gcell = .true.
  logical                      :: update_accel = .true.
  integer                      :: sizeZ, sizeY, sizeX, ib, ie, error
  real                         :: del(MDIM)
  real                         :: dvol, time, timeinfile, meandens, ekin_added, ekin_added_red
  real                         :: mass, totvol, xMomentum, yMomentum, zMomentum, xForce, yForce, zForce

  integer, parameter :: funit = 22
  character(len=80)  :: outfile = "stir.dat"

  integer, parameter :: nGlobalSum = 7                     ! Number of globally-summed quantities
  real               :: globalSumQuantities (0:nGlobalSum) ! Global summed quantities
  real               :: localSumQuantities (0:nGlobalSum)  ! Global summed quantities

  real, DIMENSION(:,:,:,:), POINTER :: solnData

#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC) :: ke, ke1, ekin_added_block
  real, dimension(GRID_IHI_GC) :: iCoord
  real, dimension(GRID_JHI_GC) :: jCoord
  real, dimension(GRID_KHI_GC) :: kCoord
#else
  real,allocatable, dimension(:) :: ke, ke1, ekin_added_block
  real,allocatable, dimension(:) :: iCoord, jCoord, kCoord
  integer :: istat
#endif

  if (firstCall) then
        call Driver_getSimTime(time)
        st_lastTimeUpdatedAccel = time-0.9999*st_dtUpdateAccel ! first time stirring
        if (st_globalMe == MASTER_PE) then
          open(funit, file=trim(outfile), position='APPEND')
          write(funit,'(10(1X,A16))') '[00]time', '[01]dt', '[02]d(Ekin)', '[03]d(Ekin)/dt', &
                                      '[04]xForce', '[05]yForce', '[06]zForce', &
                                      '[07]xMomentum', '[08]yMomentum', '[09]zMomentum'
          close(funit)
        endif
        firstCall = .false.
  endif

  ! if not using stirring, return
  if (.not.st_useStir) return

  call Timers_start("stir_timer")

  ! check if we need to update acceleration field
  update_accel = .false.
  call Driver_getSimTime(time)
  if (time .ge. st_lastTimeUpdatedAccel + st_dtUpdateAccel) then
     call st_readStirringDataFromFile(st_infilename, time, timeinfile)
     st_lastTimeUpdatedAccel = timeinfile
     update_accel = .true.
  endif

  globalSumQuantities = 0.0
  localSumQuantities  = 0.0

! use this preprocessor switch to make the Ornstein-Uhlenbeck field be treated
! as a force field (#define TREAT_AS_FORCE)
! or as an acceleration field, the default (#undef TREAT_AS_FORCE)
! Note that treating it as a force is impossible for supersonic turbulence
! (haven't tried for subsonic, but guess it will be a similar problem for finite Mach > 0),
! because it leads to local runaway accelerations in cells with low density
! and thus to diverging velocities in such cells.
#undef TREAT_AS_FORCE

#define CORRECT_BULK_MOTION
#ifdef CORRECT_BULK_MOTION

  ! sum quantities over list of blocks
  do BlockID = 1, blockCount

     !get the index limits of the block
     call Grid_getBlkIndexLimits(blockList(BlockID), blkLimits, blkLimitsGC)

     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockList(BlockID), solnData)

     !getting the dx's
     call Grid_getDeltas(blocklist(BlockID), del)

#if NDIM == 1
     dvol = del(IAXIS)
#endif
#if NDIM == 2
     dvol = del(IAXIS) * del(JAXIS)
#endif
#if NDIM == 3
     dvol = del(IAXIS) * del(JAXIS) * del(KAXIS)
#endif

     sizeX=blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
     sizeY=blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
     sizeZ=blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

#ifndef FIXEDBLOCKSIZE
     allocate(iCoord(sizeX),stat=istat)
     if (istat .ne. 0) call Driver_abortFlash("could not allocate iCoord in Stir.F90")
     allocate(jCoord(sizeY),stat=istat)
     if (istat .ne. 0) call Driver_abortFlash("could not allocate jCoord in Stir.F90")
     allocate(kCoord(sizeZ),stat=istat)
     if (istat .ne. 0) call Driver_abortFlash("could not allocate kCoord in Stir.F90")
#endif
     ! x coordinates
     call Grid_getCellCoords(IAXIS,blockList(blockID),CENTER,gcell,iCoord,sizeX)
#if NDIM > 1
     ! y coordinates
     call Grid_getCellCoords(JAXIS,blockList(blockID),CENTER,gcell,jCoord,sizeY)
#if NDIM > 2
     ! z coordinates
     call Grid_getCellCoords(KAXIS,blockList(blockID),CENTER,gcell,kCoord,sizeZ)
#endif
#endif
     ! update forcing pattern, otherwise use previous forcing pattern
     if (update_accel) call st_calcAccel(blockList(blockID),blkLimitsGC,iCoord,jCoord,kCoord)

     ! Sum contributions from the indicated blkLimits of cells.
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

              ! volume
              localSumQuantities(0) = localSumQuantities(0) + dvol

              ! mass
#ifdef DENS_VAR
              localSumQuantities(1) = localSumQuantities(1) + solnData(DENS_VAR,i,j,k)*dvol

              ! momentum
#ifdef VELX_VAR
              localSumQuantities(2) = localSumQuantities(2) + solnData(DENS_VAR,i,j,k)*solnData(VELX_VAR,i,j,k)*dvol
#endif
#ifdef VELY_VAR
              localSumQuantities(3) = localSumQuantities(3) + solnData(DENS_VAR,i,j,k)*solnData(VELY_VAR,i,j,k)*dvol
#endif
#ifdef VELZ_VAR
              localSumQuantities(4) = localSumQuantities(4) + solnData(DENS_VAR,i,j,k)*solnData(VELZ_VAR,i,j,k)*dvol
#endif
              ! driving force
#ifdef ACCX_VAR
#ifdef TREAT_AS_FORCE
              localSumQuantities(5) = localSumQuantities(5) + solnData(ACCX_VAR,i,j,k)*dvol
#else
              localSumQuantities(5) = localSumQuantities(5) + solnData(DENS_VAR,i,j,k)*solnData(ACCX_VAR,i,j,k)*dvol
#endif
#endif
#ifdef ACCY_VAR
#ifdef TREAT_AS_FORCE
              localSumQuantities(6) = localSumQuantities(6) + solnData(ACCY_VAR,i,j,k)*dvol
#else
              localSumQuantities(6) = localSumQuantities(6) + solnData(DENS_VAR,i,j,k)*solnData(ACCY_VAR,i,j,k)*dvol
#endif
#endif
#ifdef ACCZ_VAR
#ifdef TREAT_AS_FORCE
              localSumQuantities(7) = localSumQuantities(7) + solnData(ACCZ_VAR,i,j,k)*dvol
#else
              localSumQuantities(7) = localSumQuantities(7) + solnData(DENS_VAR,i,j,k)*solnData(ACCZ_VAR,i,j,k)*dvol
#endif
#endif
#endif
! ifdef DENS_VAR
           enddo
        enddo
     enddo

     call Grid_releaseBlkPtr(blockList(BlockID),solnData)
#ifndef FIXEDBLOCKSIZE
     deallocate(iCoord)
     deallocate(jCoord)
     deallocate(kCoord)
#endif

  enddo

  ! now communicate all global summed quantities to all processors
  call MPI_AllReduce (localSumQuantities, globalSumQuantities, nGlobalSum+1, &
                      FLASH_REAL, MPI_Sum, st_meshComm, error)

  totvol    = globalSumQuantities (0)
  mass      = globalSumQuantities (1)
  xMomentum = globalSumQuantities (2)
  yMomentum = globalSumQuantities (3)
  zMomentum = globalSumQuantities (4)
  xForce    = globalSumQuantities (5)
  yForce    = globalSumQuantities (6)
  zForce    = globalSumQuantities (7)

  meandens = mass / totvol

#endif
! ifdef CORRECT_BULK_MOTION

  ! set to zero for adding processor contributions below
  ekin_added = 0.

  ! Loop over local blocks
  do blockID = 1, blockCount

     ! getting the dx's
     call Grid_getDeltas(blocklist(BlockID), del)

#if NDIM == 1
     dvol = del(IAXIS)
#endif
#if NDIM == 2
     dvol = del(IAXIS) * del(JAXIS)
#endif
#if NDIM == 3
     dvol = del(IAXIS) * del(JAXIS) * del(KAXIS)
#endif

     ! Get cell coordinates for this block
     call Grid_getBlkIndexLimits(blockList(blockID),blkLimits,blkLimitsGC)

     sizeX=blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
     sizeY=blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
     sizeZ=blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

#ifndef FIXEDBLOCKSIZE
     allocate(ke(sizeX),stat=istat)
     if (istat .ne. 0) call Driver_abortFlash("could not allocate ke in Stir.F90")
     allocate(ke1(sizeX),stat=istat)
     if (istat .ne. 0) call Driver_abortFlash("could not allocate ke1 in Stir.F90")
     allocate(ekin_added_block(sizeX),stat=istat)
     if (istat .ne. 0) call Driver_abortFlash("could not allocate ekin_added_block in Stir.F90")
     allocate(iCoord(sizeX),stat=istat)
     if (istat .ne. 0) call Driver_abortFlash("could not allocate iCoord in Stir.F90")
     allocate(jCoord(sizeY),stat=istat)
     if (istat .ne. 0) call Driver_abortFlash("could not allocate jCoord in Stir.F90")
     allocate(kCoord(sizeZ),stat=istat)
     if (istat .ne. 0) call Driver_abortFlash("could not allocate kCoord in Stir.F90")
#endif

     ! set to zero for adding block and cell contributions below
     ekin_added_block(:) = 0.

! update acceleration ?
#ifndef CORRECT_BULK_MOTION
     ! x coordinates
     call Grid_getCellCoords(IAXIS,blockList(blockID),CENTER,gcell,iCoord,sizeX)
#if NDIM > 1
     ! y coordinates
     call Grid_getCellCoords(JAXIS,blockList(blockID),CENTER,gcell,jCoord,sizeY)
#if NDIM > 2
     ! z coordinates
     call Grid_getCellCoords(KAXIS,blockList(blockID),CENTER,gcell,kCoord,sizeZ)
#endif
#endif
        ! update forcing pattern, otherwise use previous forcing pattern
        if (update_accel) call st_calcAccel(blockList(blockID),blkLimitsGC,iCoord,jCoord,kCoord)
#endif
! ifndef CORRECT_BULK_MOTION

     call Grid_getBlkPtr(blockList(blockID),solnData,CENTER)
     ib = blkLimits(LOW, IAXIS)
     ie = blkLimits(HIGH, IAXIS)
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)

        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)

#if NDIM == 1
           ke(ib:ie) = 0.5*(solnData(VELX_VAR,ib:ie,j,k)**2)
#endif
#if NDIM == 2
           ke(ib:ie) = 0.5*(solnData(VELX_VAR,ib:ie,j,k)**2+solnData(VELY_VAR,ib:ie,j,k)**2)
#endif
#if NDIM == 3
           ke(ib:ie) = 0.5*(solnData(VELX_VAR,ib:ie,j,k)**2+ &
                            solnData(VELY_VAR,ib:ie,j,k)**2+ &
                            solnData(VELZ_VAR,ib:ie,j,k)**2)
#endif

#ifdef CORRECT_BULK_MOTION
#ifdef TREAT_AS_FORCE
           solnData(VELX_VAR,ib:ie,j,k) = solnData(VELX_VAR,ib:ie,j,k) + &
                (solnData(ACCX_VAR,ib:ie,j,k)/solnData(DENS_VAR,ib:ie,j,k)*meandens-xForce/totvol)*dt - &
                xMomentum/mass
#else
           solnData(VELX_VAR,ib:ie,j,k) = solnData(VELX_VAR,ib:ie,j,k) + &
                (solnData(ACCX_VAR,ib:ie,j,k)-xForce/mass)*dt - xMomentum/mass
#endif
#if NDIM > 1
#ifdef TREAT_AS_FORCE
           solnData(VELY_VAR,ib:ie,j,k) = solnData(VELY_VAR,ib:ie,j,k) + &
                (solnData(ACCY_VAR,ib:ie,j,k)/solnData(DENS_VAR,ib:ie,j,k)*meandens-yForce/totvol)*dt - &
                yMomentum/mass
#else
           solnData(VELY_VAR,ib:ie,j,k) = solnData(VELY_VAR,ib:ie,j,k) + &
                (solnData(ACCY_VAR,ib:ie,j,k)-yForce/mass)*dt - yMomentum/mass
#endif
#if NDIM > 2
#ifdef TREAT_AS_FORCE
           solnData(VELZ_VAR,ib:ie,j,k) = solnData(VELZ_VAR,ib:ie,j,k) + &
                (solnData(ACCZ_VAR,ib:ie,j,k)/solnData(DENS_VAR,ib:ie,j,k)*meandens-zForce/totvol)*dt - &
                zMomentum/mass
#else
           solnData(VELZ_VAR,ib:ie,j,k) = solnData(VELZ_VAR,ib:ie,j,k) + &
                (solnData(ACCZ_VAR,ib:ie,j,k)-zForce/mass)*dt - zMomentum/mass
#endif
#endif
#endif
#endif
! ifdef CORRECT_BULK_MOTION
#ifndef CORRECT_BULK_MOTION
           solnData(VELX_VAR,ib:ie,j,k) = solnData(VELX_VAR,ib:ie,j,k) + &
                solnData(ACCX_VAR,ib:ie,j,k)*dt
#if NDIM > 1
           solnData(VELY_VAR,ib:ie,j,k) = solnData(VELY_VAR,ib:ie,j,k) + &
                solnData(ACCY_VAR,ib:ie,j,k)*dt
#if NDIM > 2
           solnData(VELZ_VAR,ib:ie,j,k) = solnData(VELZ_VAR,ib:ie,j,k) + &
                solnData(ACCZ_VAR,ib:ie,j,k)*dt
#endif
#endif
#endif
! ifdef CORRECT_BULK_MOTION

#if NDIM == 1
           ke1(ib:ie) = 0.5*(solnData(VELX_VAR,ib:ie,j,k)**2)
#endif
#if NDIM == 2
           ke1(ib:ie) = 0.5*(solnData(VELX_VAR,ib:ie,j,k)**2+solnData(VELY_VAR,ib:ie,j,k)**2)
#endif
#if NDIM == 3
           ke1(ib:ie) = 0.5*(solnData(VELX_VAR,ib:ie,j,k)**2+ &
                             solnData(VELY_VAR,ib:ie,j,k)**2+ &
                             solnData(VELZ_VAR,ib:ie,j,k)**2)
#endif
           ! update the total energy
           solnData(ENER_VAR,ib:ie,j,k) = solnData(ENER_VAR,ib:ie,j,k) + (ke1(ib:ie)-ke(ib:ie))

           ! compute injected kinetic energy and add it to the sum
           ekin_added_block(ib:ie) = ekin_added_block(ib:ie) + &
                                      (ke1(ib:ie)-ke(ib:ie))*solnData(DENS_VAR,ib:ie,j,k)*dvol

        enddo
     enddo

     ! sum up the ekin_added block contributions
     do i = ib, ie
       ekin_added = ekin_added + ekin_added_block(i)
     enddo

     call Grid_releaseBlkPtr(blockList(blockID),solnData)
#ifndef FIXEDBLOCKSIZE
     deallocate(ke)
     deallocate(ke1)
     deallocate(ekin_added_block)
     deallocate(iCoord)
     deallocate(jCoord)
     deallocate(kCoord)
#endif
  enddo ! loop over blocks

  ! sum up injected kinetic energy contributions from all blocks and processors
  ekin_added_red = 0.
  call MPI_AllReduce (ekin_added, ekin_added_red, 1, &
       FLASH_REAL, MPI_Sum, st_meshComm, error)
  ekin_added = ekin_added_red

  ! write time evolution of ekin_added to file
  if (st_globalMe == MASTER_PE) then
    open(funit, file=trim(outfile), position='APPEND')
    write(funit,'(10(1X,ES16.9))') time, dt, ekin_added, ekin_added/dt, &
                                   xForce, yForce, zForce, xMomentum, yMomentum, zMomentum
    close(funit)
  endif

  call Timers_stop ("stir_timer")

  return

end subroutine Stir
