!!****if* source/Simulation/SimulationMain/unitTest/Grid/GCellValgrindError/Grid_unitTest
!!
!! NAME
!!
!!  Grid_unitTest
!!
!! SYNOPSIS
!!
!!  Grid_unitTest(integer(in):: fileUnit,
!!                logical(inout)::perfect)
!!
!! DESCRIPTION
!!
!!  This routine tests the internal operation of meshes
!!  Specifically it tests guardcell filling with the routine
!!  Grid_fillGuardCells.  If all tests pass then a file
!!  usually named unitTest_000x is written with the line, "All
!!  tests conformed to expected results".  If the tests failed,
!!  various error messages may be written to the file.
!!  
!!
!! ARGUMENTS
!!
!!  fileUnit : logical unit number for file in which to write error messages
!!  perfect : indicates if all tests passed or not
!!
!!***

subroutine Grid_unitTest(fileUnit,perfect)

  use Grid_interface, ONLY : Grid_fillGuardCells  
  implicit none

#include "Flash.h"
#include "constants.h"
  
  integer, intent(in)           :: fileUnit ! Output to file
  logical, intent(inout)        :: perfect  ! Flag to indicate errors

  call local_InitInternalCells()
  call local_clearGuardcells()
  call local_checkInternalCells()

  call Grid_fillGuardCells(CENTER,ALLDIR)

  call local_checkInternalCells()
  call local_checkAllCells()

  perfect=.true.

end subroutine Grid_unitTest


#define GARBAGE_VALUE -999.0

subroutine local_InitAllCells
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkPtr, &
       Grid_releaseBlkPtr
  use Grid_data, ONLY : gr_meshMe
  implicit none
  real, pointer, dimension(:,:,:,:) :: solnData
  integer, dimension(MAXBLOCKS) :: blkList
  integer :: blkCount, blk

  call Grid_getListOfBlocks(LEAF, blkList, blkCount)
  do blk = 1, blkCount
     call Grid_getBlkPtr(blkList(blk), solnData)
     solnData(:,:,:,:) = gr_meshMe
     call Grid_releaseBlkPtr(blkList(blk), solnData)
  end do
end subroutine local_InitAllCells


subroutine local_InitInternalCells
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkIndexLimits, &
       Grid_getBlkPtr, Grid_releaseBlkPtr
  use Grid_data, ONLY : gr_meshMe
  implicit none
  real, pointer, dimension(:,:,:,:) :: solnData
  integer, dimension(MAXBLOCKS) :: blkList
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer :: blkCount, blk

  call Grid_getListOfBlocks(LEAF, blkList, blkCount)
  do blk = 1, blkCount
     call Grid_getBlkIndexLimits(blkList(blk), blkLimits, blkLimitsGC)
     call Grid_getBlkPtr(blkList(blk), solnData)
     solnData(:, &
          blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = gr_meshMe
     call Grid_releaseBlkPtr(blkList(blk), solnData)
  end do
end subroutine local_InitInternalCells


subroutine local_clearGuardCells
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkIndexLimits, &
       Grid_getBlkPtr, Grid_releaseBlkPtr
  implicit none
  real, allocatable, dimension(:) :: garbage
  real, pointer, dimension(:,:,:,:) :: solnData
  integer, dimension(MAXBLOCKS) :: blkList
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer :: blkCount, blk, i, j, k

  allocate(garbage(1))
  !garbage(1) = GARBAGE_VALUE !Comment out to initialize with uninitialized heap data - useful for valgrind.

  call Grid_getListOfBlocks(LEAF, blkList, blkCount)
  do blk = 1, blkCount
     call Grid_getBlkIndexLimits(blkList(blk), blkLimits, blkLimitsGC)
     call Grid_getBlkPtr(blkList(blk), solnData)
     do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
        do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
           do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)

              if ( k < blkLimits(LOW,KAXIS) .or. k > blkLimits(HIGH,KAXIS) .or. &
                   j < blkLimits(LOW,JAXIS) .or. j > blkLimits(HIGH,JAXIS) .or. &
                   i < blkLimits(LOW,IAXIS) .or. i > blkLimits(HIGH,IAXIS) ) then
                 solnData(:,i,j,k) = garbage(1)
              end if

           end do
        end do
     end do
     call Grid_releaseBlkPtr(blkList(blk), solnData)
  end do
  deallocate(garbage)
end subroutine local_clearGuardCells


subroutine local_checkAllCells
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkPtr, &
       Grid_releaseBlkPtr
  implicit none
  real, pointer, dimension(:,:,:,:) :: solnData
  integer, dimension(MAXBLOCKS) :: blkList
  integer :: blkCount, blk
  character(len=*), parameter :: errMsg = &
       "*** Garbage in unk (test included guard cells) ***"
  logical :: isError

  isError = .false.
  call Grid_getListOfBlocks(LEAF, blkList, blkCount)
  do blk = 1, blkCount
     call Grid_getBlkPtr(blkList(blk), solnData)
     if (any(solnData(:,:,:,:) == GARBAGE_VALUE)) then
        isError = .true.
     end if
     call Grid_releaseBlkPtr(blkList(blk), solnData)
     if (isError) exit
  end do

  if (isError) then
     print *, errMsg
     call Driver_abortFlash(errMsg)
  end if
end subroutine local_checkAllCells


subroutine local_checkInternalCells
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkIndexLimits, &
       Grid_getBlkPtr, Grid_releaseBlkPtr
  implicit none
  real, pointer, dimension(:,:,:,:) :: solnData
  integer, dimension(MAXBLOCKS) :: blkList
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer :: blkCount, blk
  character(len=*), parameter :: errMsg = &
       "*** Garbage in unk (test excluded guard cells) ***"
  logical :: isError

  isError = .false.
  call Grid_getListOfBlocks(LEAF, blkList, blkCount)
  do blk = 1, blkCount
     call Grid_getBlkIndexLimits(blkList(blk), blkLimits, blkLimitsGC)
     call Grid_getBlkPtr(blkList(blk), solnData)
     if (any( solnData(:, &
          blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) == GARBAGE_VALUE )) then
        isError = .true.
     end if
     call Grid_releaseBlkPtr(blkList(blk), solnData)
     if (isError) exit
  end do

  if (isError) then
     print *, errMsg
     call Driver_abortFlash(errMsg)
  end if
end subroutine local_checkInternalCells
