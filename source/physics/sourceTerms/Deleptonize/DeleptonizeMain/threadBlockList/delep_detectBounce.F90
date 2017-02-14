!!****if* source/physics/sourceTerms/Deleptonize/DeleptonizeMain/threadBlockList/delep_detectBounce
!!
!! NAME
!!  
!!  delep_detectBounce 
!!
!!
!! SYNOPSIS
!! 
!!  call delep_detectBounce (integer(IN) :: blockCount,
!!             integer(IN) :: blockList(blockCount),
!!             real(IN)    :: dt,
!!             real(IN)    :: time)
!!
!!  
!!  
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!  blockCount : number of blocks to operate on
!!  blockList  : list of blocks to operate on
!!  dt         : current timestep
!!  time       : current time
!!
!!***

subroutine delep_detectBounce (blockCount,blockList,dt,time)
!
!==============================================================================
!
#include "Flash.h"
#include "constants.h"

  use Deleptonize_data, ONLY : delep_postBounce, delep_meshComm, delep_meshMe, &
       delep_bounceTime, delep_maxDens, delep_centralDens, delep_anteSonic, &
       delep_centralEntr
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_getBlkIndexLimits, Grid_getCellCoords
  use Logfile_interface, ONLY : Logfile_stampMessage
  !$ use omp_lib 
  use Timers_interface, ONLY : Timers_start, Timers_stop

  implicit none
  include "Flash_mpi.h"
  
  integer,intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN)::blockList
  real,intent(IN) :: dt,time
  
  real, pointer, dimension(:,:,:,:) :: solnData
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC

  integer,dimension(MDIM)  :: dimSize
  real,allocatable, dimension(:) :: xCenter, yCenter, zCenter

  integer :: blockID
  integer :: i,j,k,n
  integer :: ierr
  logical :: threadBlockList
  logical :: gcell = .true.

  character(len=MAX_STRING_LENGTH)  :: message

  real, dimension(MAXBLOCKS) :: blkMaxDens, blkMaxEntr, blkAnteS, blkMinEntr
  real :: localMaxDens, globalMaxDens, bounceDens, globalMaxAnteS
  real :: localMaxEntr, globalMaxEntr, globalMinEntr
  real :: localMaxAnteS, localMinEntr
  real :: radius
  real, dimension(3) :: localMax, globalMax

  real, allocatable, dimension(:,:,:) :: sndSpd, velEsc, factor
  integer :: imin,imax,jmin,jmax,kmin,kmax

#ifdef ST_THREAD_BLOCK_LIST
  threadBlockList = .true.

#ifdef ST_THREAD_WITHIN_BLOCK
  call Driver_abortFlash("Cannot include both threading strategies")
#endif

#else
  threadBlockList = .false.
#endif


!!$  if (delep_postBounce) return !bounce already detected

  localMaxAnteS = 0.

  !$omp parallel if(threadBlockList) &
  !$omp default(none) &
  !$omp private(n,i,j,k,blockID,blkLimits,blkLimitsGC,solnData,dimsize,xCenter,yCenter,zCenter,&
  !$omp radius,imin,imax,jmin,jmax,kmin,kmax,sndSpd,velEsc,factor) &
  !$omp shared(blockCount,blockList,localMaxDens,blkMaxDens,blkMaxEntr,blkMinEntr,blkAnteS,localMaxEntr,delep_postBounce,gcell)

  blkMaxDens = 0.
  blkAnteS = 0.
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

     allocate(sndSpd(imax-imin+1,jmax-jmin+1,kmax-kmin+1))
     allocate(velEsc(imax-imin+1,jmax-jmin+1,kmax-kmin+1))
     allocate(factor(imax-imin+1,jmax-jmin+1,kmax-kmin+1))

     blkMaxDens(n) = maxval(solnData(DENS_VAR,imin:imax,&
                               jmin:jmax,&
                               kmin:kmax))

!!     localMaxDens = max(localMaxDens, blkMaxDens)
     if (delep_postBounce) then

        factor = min(1.,max(0.,(solnData(DENS_VAR,imin:imax,jmin:jmax,kmin:kmax) - 1.0e7)))
        
        sndSpd = sqrt(solnData(GAMC_VAR,imin:imax,jmin:jmax,kmin:kmax)*solnData(PRES_VAR,imin:imax,jmin:jmax,kmin:kmax) &
                      /solnData(DENS_VAR,imin:imax,jmin:jmax,kmin:kmax)) * factor
        velEsc = sqrt(-2.*solnData(GPOT_VAR,imin:imax,jmin:jmax,kmin:kmax))

        blkAnteS(n) = maxval( (sndSpd/velEsc)**2 / solnData(GAMC_VAR,imin:imax,jmin:jmax,kmin:kmax) )

     else
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

!!        localMaxEntr = max(localMaxEntr, blkMaxEntr)

     deallocate(xCenter)
#if NDIM >1
     deallocate(yCenter)
#if NDIM == 3
     deallocate(zCenter)
#endif
#endif

     endif
     call Grid_releaseBlkPtr(blockID,solndata)
     deallocate(velEsc)
     deallocate(sndSpd)
     deallocate(factor)
  enddo
  !$omp enddo
  !$omp end parallel

  localMaxDens = maxval(blkMaxDens)
  localMaxEntr = maxval(blkMaxEntr)
  localMaxAnteS = maxval(blkAnteS)
  localMinEntr = minval(blkMinEntr)

  localMax(1:3) = (/localMaxDens,localMaxEntr,localMaxAnteS/)

  call Timers_start("allReduce")
  call MPI_Allreduce(localMax, globalMax, 3, FLASH_REAL, MPI_MAX, &
                     delep_meshComm, ierr)
  call MPI_Allreduce(localMinEntr, globalMinEntr, 1, FLASH_REAL, MPI_MIN, &
                     delep_meshComm, ierr)
  call Timers_stop("allReduce")

  delep_centralDens = globalMax(1)
  globalMaxEntr = globalMax(2)
  delep_anteSonic = globalMax(3)
  delep_centralEntr = globalMinEntr

!!$  call MPI_Allreduce(localMaxDens, globalMaxDens, 1, FLASH_REAL, MPI_MAX, &
!!$                     delep_meshComm, ierr)
!!$
!!$  delep_centralDens = globalMaxDens
!!$
!!$  call MPI_Allreduce(localMaxAnteS, globalMaxAnteS, 1, FLASH_REAL, MPI_MAX, &
!!$                     delep_meshComm, ierr)
!!$  delep_anteSonic = globalMaxAnteS

  if (.NOT. delep_postBounce) then
!!$     call MPI_Allreduce(localMaxEntr, globalMaxEntr, 1, FLASH_REAL, MPI_MAX, &
!!$                        delep_meshComm, ierr)
     if (globalMaxEntr > 3. .AND. delep_centralDens > 2.e14) then
        ! Bounce
        delep_postBounce = .TRUE.
        delep_bounceTime = time
        if (delep_meshMe == MASTER_PE) then
           print *, "Bounce!", time, delep_centralDens, globalMaxEntr
           write(message,*) "Core Bounce! Time = ", time, delep_centralDens
           call Logfile_stampMessage(message)
        endif
     endif
!!$     if (delep_maxDens > 1.e14 .AND. globalMaxDens < delep_maxDens) then
!!$        ! Bounce
!!$        delep_postBounce = .TRUE.
!!$        delep_bounceTime = time
!!$        if (delep_meshMe == MASTER_PE) then
!!$           print *, "Bounce!", time, delep_maxDens, bounceDens
!!$           write(message,*) "Core Bounce! Time = ", time, delep_maxDens
!!$           call Logfile_stampMessage(message)
!!$        endif
!!$        delep_maxDens = globalMaxDens
!!$     else
!!$        delep_maxDens = globalMaxDens
!!$     endif
  endif

!!$  call MPI_Allreduce(local_postBounce, delep_postBounce, 1, MPI_LOGICAL, &
!!$                     MPI_LOR, delep_meshComm, ierr)
!!$
!!$  delep_bounceTime = time
!!$  if (delep_postBounce .AND. delep_meshMe == MASTER_PE) then
!!$     print *, "Bounce!", time
!!$     write(message,*) "Core Bounce! Time = ", time
!!$     call Logfile_stampMessage(message)
!!$  endif

  return
end subroutine delep_detectBounce
