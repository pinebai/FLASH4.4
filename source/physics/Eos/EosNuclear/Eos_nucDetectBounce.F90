!!****if* source/physics/Eos/EosNuclear/Eos_nucDetectBounce
!!
!! NAME
!!  
!!  Eos_nucDetectBounce 
!!
!!
!! SYNOPSIS
!! 
!!  call Eos_nucDetectBounce(logical(OUT)  :: postbounce,
!!                           real, optional(OUT)  :: bouncetime,
!!                           real, optional(OUT)  :: centraldens,
!!                           real, optional(OUT)  :: centralentr)
!!
!!  
!!  
!! DESCRIPTION
!!  This routine determines if collapse has proceeded to the point of 
!!  core bounce, as determined by the maximum density.
!!
!! ARGUMENTS
!!
!!   postbounce : 
!!
!!   bouncetime : 
!!
!!   centraldens : 
!!
!!   centralentr : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

subroutine Eos_nucDetectBounce(postBounce,bounceTime,centralDens,centralEntr)
  !
  !==============================================================================
  !
#include "Flash.h"
#include "constants.h"

  use eos_nucData, ONLY : eos_postBounce, eos_bounceTime, &
       eos_centralDens, eos_centralEntr, eos_nstep, &
       eos_meshComm, eos_meshMe, eos_bounceDens, eos_shockEntr
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_getBlkIndexLimits, Grid_getCellCoords, &
       Grid_getListOfBlocks
  use Logfile_interface, ONLY : Logfile_stampMessage
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_getNStep, Driver_getSimTime
  use IO_interface, ONLY : IO_setScalar

  implicit none
  include "Flash_mpi.h"

  logical, intent(OUT) :: postBounce
  real, optional, intent(OUT) :: bounceTime, centralDens, centralEntr

  integer :: blockCount
  integer :: blockList(MAXBLOCKS)

  real, pointer, dimension(:,:,:,:) :: solnData
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC

  integer,dimension(MDIM)  :: dimSize
  real,allocatable, dimension(:) :: xCenter, yCenter, zCenter

  integer :: blockID
  integer :: i,j,k,n
  integer :: ierr
  logical :: threadBlockList
  logical :: gcell = .true.

  character(len=100)  :: message

  real, dimension(MAXBLOCKS) :: blkMaxDens, blkMaxEntr, blkMinEntr
  real :: localMaxDens
  real :: localMaxEntr, globalMaxEntr, globalMinEntr
  real :: localMinEntr
  real :: radius
  real, dimension(2) :: localMax, globalMax

  real, allocatable, dimension(:,:,:) :: sndSpd, velEsc, factor
  integer :: imin,imax,jmin,jmax,kmin,kmax

  integer :: nstep
  real :: time

#ifdef ST_THREAD_BLOCK_LIST
  threadBlockList = .true.

#ifdef ST_THREAD_WITHIN_BLOCK
  call Driver_abortFlash("Cannot include both threading strategies")
#endif

#else
  threadBlockList = .false.
#endif


  if (eos_postBounce) then
     postBounce = .TRUE.
     if (present(bounceTime)) bounceTime = eos_bounceTime
     if (present(centralDens)) centralDens = eos_centralDens
     if (present(centralEntr)) centralEntr = eos_centralEntr
     return !bounce already detected
  end if

  ! Calling routine would like us to check for bounce!
  ! First, verify that we haven't done this already for
  ! this time step.
  call Driver_getNStep(nstep)
  if (nstep == eos_nstep) then
     ! We've already checked for bounce
     postBounce = .FALSE.
     if (present(bounceTime)) bounceTime = 0.0
     if (present(centralDens)) centralDens = eos_centralDens
     if (present(centralEntr)) centralEntr = eos_centralEntr
     return
  end if

  ! We have yet to check for bounce this step
  eos_nstep = nstep

  ! Now proceed with check
  call Driver_getSimTime(time)

  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

  !$omp parallel if(threadBlockList) &
  !$omp default(none) &
  !$omp private(n,i,j,k,blockID,blkLimits,blkLimitsGC,solnData,dimsize,xCenter,yCenter,zCenter,&
  !$omp radius,imin,imax,jmin,jmax,kmin,kmax,sndSpd,velEsc,factor) &
  !$omp shared(blockCount,blockList,localMaxDens,blkMaxDens,blkMaxEntr,blkMinEntr,localMaxEntr,eos_postBounce,gcell)

  blkMaxDens = 0.
  blkMaxEntr = 0.
  blkMinEntr = 100.
  localMaxEntr = 0.
  localMaxDens = 0.

  !$omp do schedule(static) 
  do n = 1, blockCount
     blockID = blockList(n)

     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blockID,solnData)

     imin = blkLimits(LOW,IAXIS)
     imax = blkLimits(HIGH,IAXIS)
     jmin = blkLimits(LOW,JAXIS)
     jmax = blkLimits(HIGH,JAXIS)
     kmin = blkLimits(LOW,KAXIS)
     kmax = blkLimits(HIGH,KAXIS)

     blkMaxDens(n) = maxval(solnData(DENS_VAR,imin:imax,&
          jmin:jmax,&
          kmin:kmax))

     dimSize(:)=blkLimitsGC(HIGH,:)-blkLimitsGC(LOW,:)+1
     if (NDIM > 2)then
        allocate(zCenter(dimSize(KAXIS)))
        call Grid_getCellCoords(KAXIS,blockID,&
             CENTER,gcell,zCenter,dimSize(KAXIS))
     end if
     if (NDIM > 1)then
        allocate(yCenter(dimSize(JAXIS)))
        call Grid_getCellCoords(JAXIS,blockID,&
             CENTER,gcell,yCenter,dimSize(JAXIS))
     end if

     allocate(xCenter(dimSize(IAXIS)))
     call Grid_getCellCoords(IAXIS,blockID,&
          CENTER,gcell,xCenter,dimSize(IAXIS))

     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

              radius = xCenter(i)**2
#if NDIM > 1
              radius = radius + yCenter(j)**2
#if NDIM == 3
              radius = radius + zCenter(k)**2
#endif
#endif
              radius = sqrt(radius)
              if (radius < 3.e6) blkMaxEntr(n) = max(blkMaxEntr(n), solnData(ENTR_VAR,i,j,k))
              if (radius < 3.e6) blkMinEntr(n) = min(blkMinEntr(n), solnData(ENTR_VAR,i,j,k))
           enddo
        enddo
     enddo

     deallocate(xCenter)
#if NDIM >1
     deallocate(yCenter)
#if NDIM == 3
     deallocate(zCenter)
#endif
#endif

     call Grid_releaseBlkPtr(blockID,solndata)
  enddo
  !$omp enddo
  !$omp end parallel

  localMaxDens = maxval(blkMaxDens)
  localMaxEntr = maxval(blkMaxEntr)
  localMinEntr = minval(blkMinEntr)

  localMax(1:2) = (/localMaxDens,localMaxEntr/)

  call MPI_Allreduce(localMax, globalMax, 2, FLASH_REAL, MPI_MAX, &
       eos_meshComm, ierr)
  call MPI_Allreduce(localMinEntr, globalMinEntr, 1, FLASH_REAL, MPI_MIN, &
       eos_meshComm, ierr)

  eos_centralDens = globalMax(1)
  globalMaxEntr = globalMax(2)
  eos_centralEntr = globalMinEntr

  if (globalMaxEntr > eos_shockEntr .AND. eos_centralDens > eos_bounceDens) then
     ! Bounce
     eos_postBounce = .TRUE.
     eos_bounceTime = time
     call IO_setScalar("postBounce", eos_postBounce)
     call IO_setScalar("bounceTime", eos_bounceTime)
     if (eos_meshMe == MASTER_PE) then
        write(*,*) "Bounce!", time, eos_centralDens/1e14, eos_centralEntr
        write(message,*) "Bounce!", time, eos_centralDens/1e14, eos_centralEntr
        call Logfile_stampMessage(message)
     endif
  endif

  ! return the dummy variables
  postBounce = eos_postBounce
  if (present(bounceTime)) bounceTime = eos_bounceTime
  if (present(centralDens)) centralDens = eos_centralDens
  if (present(centralEntr)) centralEntr = eos_centralEntr

  return
end subroutine Eos_nucDetectBounce
