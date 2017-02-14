!!****if* source/Simulation/SimulationMain/FlameChannel/Simulation_initRestart
!!
!! NAME
!!  Simulation_initRestart
!!
!! SYNOPSIS
!!  call Simulation_initRestart()
!!
!! DESCRIPTION
!!  This is where the user should place code for a setup that needs to adjust
!!  data on a restart, particularly if grid data, grid metadata or particle 
!!  data needs to be changed on restarting.
!!
!! ARGUMENTS
!!
!!***

subroutine Simulation_initRestart()

  use Simulation_data
  use Logfile_interface, ONLY : Logfile_stampMessage
  use Grid_interface, ONLY : Grid_getListOfBlocks, &
    Grid_getBlkIndexLimits, Grid_getPointData, &
    Grid_putPointData, Grid_getCellCoords, &
    Grid_getBlkPtr, Grid_releaseBlkPtr
  use IO_data, ONLY : io_useCollectiveHDF5

  implicit none

#include "constants.h"
#include "Flash.h"
#include "mpif.h"

  character(len=MAX_STRING_LENGTH)  :: logMesg

  integer :: blockCount, maxCount
  integer :: blockList(MAXBLOCKS)
  integer :: n, bid, ierr, changed_cells, global_cells

  real, pointer, dimension(:,:,:,:) :: solnData

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM) :: cell
  integer :: i, j, k, ext_i, ext_j, ext_k
  integer :: isize, jsize, ksize

  real, allocatable, dimension(:) :: iCenter, jCenter, kCenter
  real, allocatable, dimension(:) :: iLeft, jLeft, kLeft
  real, allocatable, dimension(:) :: iRight, jRight, kRight

  real, allocatable, dimension(:,:,:) :: velx, vely, velz
  real :: vx_l, vx_r, vx_c, vy_l, vy_c, vy_r, vz_l, vz_c, vz_r
  real :: zctrVortex, vortex_stream1, vortex_stream2
  real :: y_vort1, y_vort2, z_vort
  real :: kine, eint, flam

  real, parameter :: flam_min = 1.e-6

  logical :: use_collective

  if (.not. sim_restartVortex) return

  logMesg = "[Simulation_initRestart] opening '" //&
     trim(sim_turbfield_filename) // "'"
  call Logfile_stampMessage(logMesg)

  use_collective = io_useCollectiveHDF5
  changed_cells = 0

  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

  ! call c routines to open file and setup
  call sim_turb_field_setup( sim_turbfield_bbox(IAXIS,LOW),  &
                             sim_turbfield_bbox(IAXIS,HIGH), &
                             sim_turbfield_bbox(JAXIS,LOW),  &
                             sim_turbfield_bbox(JAXIS,HIGH), &
                             sim_turbfield_bbox(KAXIS,LOW),  &
                             sim_turbfield_bbox(KAXIS,HIGH), &
                             sim_smooth_level, sim_turbfield_filename, &
                             use_collective )
  ! the setup routine might change use_collective
  if ( use_collective ) then
     call MPI_Allreduce(blockCount, maxCount, 1, MPI_Integer, &
        MPI_MAX, MPI_COMM_WORLD, ierr)
  else
     maxCount = blockCount
  endif

  do n = 1, blockCount
     bid = blockList(n)

     call Grid_getBlkPtr(bid,solnData)

     call Grid_getBlkIndexLimits(bid, blkLimits, blkLimitsGC)

     isize = blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1
     allocate(iCenter(isize))
!     allocate(iLeft(isize))
!     allocate(iRight(isize))
     jsize = blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1
     allocate(jCenter(jsize))
!     allocate(jLeft(jsize))
!     allocate(jRight(jsize))
     ksize = blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1
     allocate(kCenter(ksize))
!     allocate(kLeft(ksize))
!     allocate(kRight(ksize))
     call Grid_getCellCoords(IAXIS,bid,CENTER,.false.,iCenter,isize)
!     call Grid_getCellCoords(IAXIS,bid,LEFT_EDGE,.false.,iLeft,isize)
!     call Grid_getCellCoords(IAXIS,bid,RIGHT_EDGE,.false.,iRight,isize)
     call Grid_getCellCoords(JAXIS,bid,CENTER,.false.,jCenter,jsize)
!     call Grid_getCellCoords(JAXIS,bid,LEFT_EDGE,.false.,jLeft,jsize)
!     call Grid_getCellCoords(JAXIS,bid,RIGHT_EDGE,.false.,jRight,jsize)
     call Grid_getCellCoords(KAXIS,bid,CENTER,.false.,kCenter,ksize)
!     call Grid_getCellCoords(KAXIS,bid,LEFT_EDGE,.false.,kLeft,ksize)
!     call Grid_getCellCoords(KAXIS,bid,RIGHT_EDGE,.false.,kRight,ksize)

     allocate(velx(isize,jsize,ksize))
     allocate(vely(isize,jsize,ksize))
     allocate(velz(isize,jsize,ksize))

     call sim_turb_field_get_vel(isize, jsize, ksize, &
        velx, vely, velz, iCenter, jCenter, kCenter)

     !-----------------------------------------------
     ! loop over all zones and init
     !-----------------------------------------------
     do k = 1, ksize
        ext_k = blkLimits(LOW,KAXIS) + k - 1
        do j = 1, jsize
           ext_j = blkLimits(LOW,JAXIS) + j - 1
           do i = 1, isize
              ext_i = blkLimits(LOW,IAXIS) + i - 1

              flam = solnData(FLAM_MSCALAR,ext_i,ext_j,ext_k)

              if ( flam < flam_min ) then

                 solnData(VELX_VAR,ext_i,ext_j,ext_k) = & 
                    solnData(VELX_VAR,ext_i,ext_j,ext_k) + sim_vrms*velx(i,j,k)
                 solnData(VELY_VAR,ext_i,ext_j,ext_k) = &
                    solnData(VELY_VAR,ext_i,ext_j,ext_k) + sim_vrms*vely(i,j,k)
                 solnData(VELZ_VAR,ext_i,ext_j,ext_k) = &
                    solnData(VELZ_VAR,ext_i,ext_j,ext_k) + sim_vrms*velz(i,j,k)

                 kine = 0.5 * ( solnData(VELX_VAR,ext_i,ext_j,ext_k)**2 + &
                                solnData(VELY_VAR,ext_i,ext_j,ext_k)**2 + &
                                solnData(VELZ_VAR,ext_i,ext_j,ext_k)**2 )

                 solnData(ENER_VAR,ext_i,ext_j,ext_k) = &
                    solnData(EINT_VAR,ext_i,ext_j,ext_k) + kine

                 changed_cells = changed_cells + 1

              endif
           enddo
        enddo
     enddo
   
     call Grid_releaseBlkPtr(bid, solnData)

     deallocate(velx)
     deallocate(vely)
     deallocate(velz)

     deallocate(iCenter)
!     deallocate(iLeft)
!     deallocate(iRight)
     deallocate(jCenter)
!     deallocate(jLeft)
!     deallocate(jRight)
     deallocate(kCenter)
!     deallocate(kLeft)
!     deallocate(kRight)

  enddo

  do n = blockCount+1, maxCount
     call sim_turb_field_null_read()
  enddo

  call sim_turb_field_teardown

  call MPI_Reduce(changed_cells, global_cells, 1, MPI_INTEGER, MPI_SUM, &
     MASTER_PE, MPI_COMM_WORLD, ierr)
  write (unit=logMesg,fmt='(i10)') global_cells
  logMesg = "[Simulation_initRestart] Modified velocity in " // &
       trim(logMesg) // " cells"
  call Logfile_stampMessage(logMesg)

  return

end subroutine Simulation_initRestart
