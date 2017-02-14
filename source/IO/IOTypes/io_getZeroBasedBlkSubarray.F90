!!****if* source/IO/IOTypes/io_getZeroBasedBlkSubarray
!!
!! NAME
!!  io_getZeroBasedBlkSubarray
!!
!! SYNOPSIS
!!
!!  io_getZeroBasedBlkSubarray(integer(IN) :: gridDataStruct,
!!                             integer(OUT) :: blockInnerSize(MDIM)
!!                             integer(OUT) :: blockOuterSize(MDIM)
!!                             integer(OUT) :: blockInnerOffset(MDIM))
!!
!!
!! DESCRIPTION
!!
!! This subroutine returns three arrays that are used to describe
!! the internal region of a block for the specified grid data structure.
!!
!!
!! ARGUMENTS
!!
!! gridDataStruct: The grid data structure, e.g. UNK, FACEX
!! blockInnerSize: The size of the internal region
!! blockOuterSize: The size of the internal region + guardcells
!! blockInnerOffset: The first cell of the internal region
!!
!!***


subroutine io_getZeroBasedBlkSubarray(gridDataStruct, blockInnerSize, &
     blockOuterSize, blockInnerOffset)

#include "constants.h"
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getBlkIndexLimits

  implicit none
  integer, intent(IN) :: gridDataStruct
  integer, dimension(MDIM), intent(OUT) :: blockInnerSize, blockOuterSize, &
       blockInnerOffset
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
        
  !We assume that each block belonging to a single process has the same size.
  !However, we allow different processes to have different block sizes.
  call Grid_getBlkIndexLimits(1, blkLimits, blkLimitsGC, gridDataStruct)
        
  blockInnerSize(1:MDIM) = &
       blkLimits(HIGH,1:MDIM) - blkLimits(LOW,1:MDIM) + 1
  blockOuterSize(1:MDIM) = &
       blkLimitsGC(HIGH,1:MDIM) - blkLimitsGC(LOW,1:MDIM) + 1
  blockInnerOffset(1:MDIM) = blkLimits(LOW,1:MDIM) - 1 !Make offset zero-based.
  
  !SCRATCH is handled specially.  The internal data is (nxb+1,nyb+1,nzb+1) 
  !but we are only interested in (nxb,nyb,nzb) at the current time.
  !NOTE: This must be consistent with io_ncmpi_write_grid_header.c.
  if (gridDataStruct == SCRATCH) then
     blockInnerSize(1:MDIM) = blockInnerSize(1:MDIM) - 1
  end if

end subroutine io_getZeroBasedBlkSubarray
