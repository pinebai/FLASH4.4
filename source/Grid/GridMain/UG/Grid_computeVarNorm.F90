!!****if* source/Grid/GridMain/UG/Grid_computeVarNorm
!!
!! NAME
!!  Grid_computeVarNorm
!!
!! SYNOPSIS
!!  
!!  call Grid_computeVarNorm(integer(in)  :: level,
!!                           integer(in)  :: normType,
!!                           integer(in)  :: ivar,
!!                           real(out)    :: norm,
!!                           integer(in)  :: leafOnly)
!!
!! DESCRIPTION
!!
!!  Computes the L1 or L2 norm of the variable specified by ivar.  This
!!  can be done per-level, or on leaf or all nodes.  For multigrid, the
!!  L2 norm is used for convergence, but the L1 norm is incredibly useful
!!  for debugging purposes.
!!
!! ARGUMENTS
!!
!!  level     - If the norm is restricted to a given level; 0 is all
!!  normType - p in the Lp norm where choices of p are 1 or 2
!!  ivar      - the grid variable being normed; -1 for work
!!  norm      - the variable with which to return the norm
!!  leafOnly - if this isn't 0, compute the norm only on leaf nodes
!!
!! RESULT
!!
!!  The norm of ivar is in norm.
!!
!! EXAMPLE
!!  
!!  gr_restrictTree()
!!  do i = 1, lrefine_max
!!    call Grid_computeVarNorm(i, 1, pdens, norm(i), 0)
!!  enddo
!!  do i = 1, lrefine_max
!!    if (norm(0) - norm(i) > 0.0000001) then
!!    driver_abortFlash("restriction is highly nonconservatory!")
!!    endif
!!  enddo
!!
!!***

!!REORDER(5): unk

subroutine Grid_computeVarNorm (level, normType, ivar, norm, leafOnly)


#include "constants.h"

  use physicaldata, ONLY : unk
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkData
  use Driver_interface, ONLY : Driver_abortFlash
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_data, ONLY : gr_meshComm

  implicit none

  include "Flash_mpi.h"

  integer, intent(IN)  :: normType, level, ivar, leafOnly
  real, intent(OUT)    :: norm
  
  integer :: i, j, k, ierr
  real    :: lvol, lsum, bsum, sum
  real    :: cvol
  logical :: include_in_sum
  integer :: totalblockshere
  integer, dimension(2,MDIM)   :: blkLimitsGC, blkLimits
  real, allocatable :: cellVolumes(:,:,:)


!===============================================================================

  call Timers_start("Grid_computeVarNorm")

  lvol = 0.
  lsum = 0.
  totalblockshere = 0

  if (normType /= 1 .and. normType /= 2) then
    call Driver_abortFlash('only L1 and L2 norms supported in Grid_computeVarNorm!')
  endif

  call Grid_getBlkIndexLimits(1,blkLimits,blkLimitsGC)
  allocate(cellVolumes(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
       blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
       blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)))

  include_in_sum = ((level == 1) .or. (level == 0))
  ! leafOnly must be ignored for UG
  if (include_in_sum) then
     totalblockshere = totalblockshere + 1
     call Grid_getBlkData(1, CELL_VOLUME, 0, EXTERIOR, &
          (/blkLimits(LOW,IAXIS),blkLimits(LOW,JAXIS),blkLimits(LOW,KAXIS)/), &
          cellVolumes, &
          (/blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1, &
            blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1, &
            blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1/) )
     bsum = 0.
     if (ivar >= 0) then
        do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)  ! working on interior only
           do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
                 cvol = cellVolumes(i,j,k)
                 bsum = bsum + abs(unk(ivar,i,j,k,1))**normType * cvol
                 lvol = lvol + cvol
              enddo
           enddo
        enddo
!!     else
        ! DEV: Issue warning?
     endif
     lsum = lsum + bsum
  endif

  deallocate(cellVolumes)
  
  call mpi_allreduce ( lsum, sum, 1, FLASH_REAL, & 
       MPI_SUM, gr_meshComm, ierr )
  !call mpi_allreduce ( lvol, vol, 1, FLASH_REAL,
  !                     MPI_SUM, FLASH_COMM, ierr )
  if (normType == 2) then
    norm = sqrt(sum)
  else if (normType == 1) then
    norm = sum
  endif

  call Timers_stop("Grid_computeVarNorm")

  !=================================================================
  
  return
end subroutine Grid_computeVarNorm
