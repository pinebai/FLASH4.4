!!****if* source/physics/sourceTerms/Stir/StirMain/Generate/Stir
!!
!! NAME
!!
!!  Stir
!!
!! SYNOPSIS
!!
!!  Stir(integer(IN)::blockCount,
!!       integer(IN)::blockList(blockCount),
!!       real(IN) :: dt)
!!
!! DESCRIPTION
!!   Apply the stirring opperator on the list of blocks provided as input
!!
!! ARGUMENTS
!!   blockCount   : The number of blocks in the list
!!   blockList(:) : The list of blocks on which to apply the stirring operator
!!   dt           : the current timestep
!!
!! NOTES
!!   
!!
!!***

!#define REPRODUCE_SMALL_EOS_SIDE_EFFECTS

!!REORDER(4): solnData

subroutine Stir(blockCount,blockList,dt)
 
  use Stir_data, ONLY : st_nmodes, st_OUphases, &
       st_freq, st_OUphases, st_OUvar, st_decay, st_randseed, st_useStir, &
       st_meshComm, &
       st_eosMode
  use Driver_interface, ONLY : Driver_getNStep, Driver_abortFlash
  use Eos_interface, ONLY : Eos_wrapped
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, &
    Grid_getDeltas, Grid_releaseBlkPtr, Grid_getCellCoords, &
    Grid_getRowData, Grid_putRowData
  use ut_randomInterface, ONLY : ut_randomSeed
 
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer,intent(IN)                       :: blockCount
  integer,dimension(blockCount), intent(IN):: blockList
  real,intent(IN)                          :: dt


  integer                      :: blockID, i, j, k, nphases, nstep
  real,dimension(NDIM)         :: a
  integer, dimension(2,MDIM)   :: blkLimits,blkLimitsGC
  logical                      :: gcell = .true., update_accel = .true.
  integer                      :: sizeZ,sizeY,sizeX
  integer :: ib,ie
  real  :: mass  
  real  :: xMomentum, yMomentum, zMomentum
  integer                                       :: count, error


  integer, parameter ::  nGlobalSum = 4          ! Number of globally-summed quantities
  real :: globalSumQuantities (nGlobalSum) !Global summed quantities
  real :: localSumQuantities (nGlobalSum) !Global summed quantities

  real :: del(MDIM)
  integer :: lb
  real    :: dvol
  integer :: point(MDIM)
  real, DIMENSION(:,:,:,:), POINTER :: solnData

#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC) ::  ke,ke1 
  real, dimension(GRID_IHI_GC) :: iCoord
  real, dimension(GRID_JHI_GC) :: jCoord
  real, dimension(GRID_KHI_GC) :: kCoord
#else
  real,allocatable, dimension(:) :: ke,ke1 
  real,allocatable, dimension(:) :: iCoord,jCoord,kCoord
  integer :: istat
#endif

  ! If not using stirring, return
  if (.not.st_useStir) then 
     return
  endif

  call Timers_start("stir_timer")

  nphases = 6* st_nmodes

  ! Zero the mass and momentum sums
  globalSumQuantities = 0.0
  localSumQuantities = 0.0
  
  ! Only update acceleration every st_freq timesteps
  update_accel = .false.

  call Driver_getNStep (nstep)

  if ( (nstep .eq. 1) .or. (mod (nstep, st_freq) .eq. 0)  ) then
     update_accel = .true. 
  endif

  ! Sum quantities over list of blocks.  
  do lb = 1, blockCount
     !get the index limits of the block
     call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)

     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockList(lb), solnData)


     !getting the dx's
     call Grid_getDeltas(blocklist(lb), del)

     
     !! This calculation is specific to cartesian co-ordinates
     !! Eventually it should get the volumes from Grid

     if(NDIM == 1) then
        dvol = del(IAXIS)
     else if(NDIM == 2) then
        dvol = del(IAXIS) * del(JAXIS)
     else
        dvol = del(IAXIS) * del(JAXIS) * del(KAXIS)
     end if

     ! Sum contributions from the indicated blkLimits of cells.
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              
              point(IAXIS) = i
              point(JAXIS) = j
              point(KAXIS) = k

              !! This is the hook for getting cell volumes for
              !! other coordinate systems.
              !! Eventually should be replaced with another more 
              !! robust call
              !! call Grid_getSingleCellVol(blockList(i), EXTERIOR, point, dvol)
     
              ! mass   
#ifdef DENS_VAR
              localSumQuantities(1) = localSumQuantities(1) + solnData(DENS_VAR,i,j,k)*dvol 
              
#ifdef VELX_VAR      
              ! momentum
              localSumQuantities(2) = localSumQuantities(2) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELX_VAR,i,j,k)*dvol 
#endif
#ifdef VELY_VAR      
              localSumQuantities(3) = localSumQuantities(3) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELY_VAR,i,j,k)*dvol
#endif
#ifdef VELZ_VAR      
              localSumQuantities(4) = localSumQuantities(4) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELZ_VAR,i,j,k)*dvol
#endif

#endif ! ifdef DENS_VAR
           
           enddo
        enddo
     enddo
     call Grid_releaseBlkPtr(blockList(lb),solnData)
  enddo

  ! Now communicate all global summed quantities to all processors
  call MPI_AllReduce (localSumQuantities, globalSumQuantities, nGlobalSum, &
       FLASH_REAL, MPI_Sum, st_meshComm, error)
  
  mass      = globalSumQuantities (1)
  xMomentum = globalSumQuantities (2)
  yMomentum = globalSumQuantities (3)
  zMomentum = globalSumQuantities (4)

  ! Loop over local blocks
  do blockID = 1, blockCount

     ! Get cell coordinates for this block
     call Grid_getBlkIndexLimits(blockList(blockID),blkLimits,blkLimitsGC)
     sizeX=blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
     sizeY=blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
     sizeZ=blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

#ifndef FIXEDBLOCKSIZE

     allocate(ke(sizeX),stat=istat)
     if (istat .ne. 0) call Driver_abortFlash("could not allocate ke in Stir.F90")

     allocate(ke1(sizeX),stat=istat)
     if (istat .ne. 0) call Driver_abortFlash("could not allocate ke in Stir.F90")

     allocate(iCoord(sizeX),stat=istat)
     if (istat .ne. 0) call Driver_abortFlash("could not allocate iCoord in Stir.F90")

     allocate(jCoord(sizeY),stat=istat)
     if (istat .ne. 0) call Driver_abortFlash("could not allocate jCoord in Stir.F90")

     allocate(kCoord(sizeZ),stat=istat)
     if (istat .ne. 0) call Driver_abortFlash("could not allocate kCoord in Stir.F90")
#endif
     ! x coordinates
     call Grid_getCellCoords(IAXIS,blockList(blockID),&
          CENTER,gcell,iCoord,sizeX)

#if NDIM > 1
     ! y coordinates
     call Grid_getCellCoords(JAXIS,blockList(blockID),&
          CENTER,gcell,jCoord,sizeY)

#if NDIM > 2
     ! z coordinates
     call Grid_getCellCoords(KAXIS,blockList(blockID),&
          CENTER,gcell,kCoord,sizeZ)
#endif
#endif

     ! Update acceleration
     
     if (update_accel) then
        call st_calcAccel (blockList(blockID),blkLimitsGC,iCoord,jCoord,kCoord)
     endif
     call Grid_getBlkPtr(blockList(blockID),solnData,CENTER)
     ib = blkLimitsGC (LOW, IAXIS)
     ie = blkLimitsGC (HIGH, IAXIS)
     do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
        
        do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)

           ke = 0.5*(solnData(VELX_VAR,ib:ie,j,k)**2 + solnData(VELY_VAR,ib:ie,j,k)**2 &
                   + solnData(VELZ_VAR,ib:ie,j,k)**2)

           solnData(VELX_VAR,ib:ie,j,k) = solnData(VELX_VAR,ib:ie,j,k) + &
                solnData(ACCX_VAR,ib:ie,j,k) * dt - xMomentum / mass
#if NDIM > 1
           solnData(VELY_VAR,ib:ie,j,k) = solnData(VELY_VAR,ib:ie,j,k) + &
                solnData(ACCY_VAR,ib:ie,j,k) * dt - yMomentum / mass
#if NDIM > 2
           solnData(VELZ_VAR,ib:ie,j,k) = solnData(VELZ_VAR,ib:ie,j,k) + &
                solnData(ACCZ_VAR,ib:ie,j,k) * dt - zMomentum / mass
#endif
#endif
           ke1 = 0.5*(solnData(VELX_VAR,ib:ie,j,k)**2 + solnData(VELY_VAR,ib:ie,j,k)**2 &
                   + solnData(VELZ_VAR,ib:ie,j,k)**2)

           solnData(ENER_VAR,ib:ie,j,k) = solnData(ENER_VAR,ib:ie,j,k) + (ke1 - ke)
           
        enddo
     enddo
#ifdef REPRODUCE_SMALL_EOS_SIDE_EFFECTS
     call Eos_wrapped(st_eosMode,blkLimits,blockList(blockID))
#endif
     call Grid_releaseBlkPtr(blockList(blockID),solnData)
#ifndef FIXEDBLOCKSIZE
     deallocate(ke)
     deallocate(ke1)
     deallocate(iCoord)
     deallocate(jCoord)
     deallocate(kCoord)
#endif
  enddo

  if (update_accel) then 
     call st_ounoiseupdate(nphases, st_OUphases, st_OUvar, &
          st_freq * dt, st_decay)
     call st_calcPhases()
     !! Store random seed in memory for later checkpoint.
     call ut_randomSeed (ut_get = st_randseed)
  endif
  
  call Timers_stop ("stir_timer")
  
  return
end subroutine Stir
