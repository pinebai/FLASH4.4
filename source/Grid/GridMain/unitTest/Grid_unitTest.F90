!!****if* source/Grid/GridMain/unitTest/Grid_unitTest
!!
!! NAME
!!
!!  Grid_unitTest
!!
!! SYNOPSIS
!!
!!  call Grid_unitTest(integer(in)    :: fileUnit,
!!                     logical(inout) :: perfect)
!!
!! DESCRIPTION
!!
!!  This routine tests the internal operation of meshes.
!!  Specifically, it tests guard cell filling with the routine
!!  Grid_fillGuardCells.  It also tests getting and putting data into
!!  the unk data structure with the routine
!!  Grid_get(Point/Row/Plane/Blk)Data and
!!  Grid_put(Point/Row/Plane/Blk)Data.  If all tests pass then a file,
!!  usually named unitTest_000x, is written with the line
!!  "All tests conformed to expected results".  If the tests failed,
!!  various error messages may be written to the file.
!!  
!!
!! ARGUMENTS
!!
!!  fileUnit : logical unit number for file in which to write error messages
!!  perfect : indicates if all tests passed or not
!!
!! NOTES
!!  AD: Aug 2011: augmented the routine to also verify that given a position
!!      coordinate, the right cell in the domain can be identified
!!
!!***

subroutine Grid_unitTest(fileUnit,perfect)

  use Grid_interface, ONLY : Grid_fillGuardCells, &
       Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getListOfBlocks, &
       Grid_getBlkIDFromPos, Grid_getSingleCellCoords, Grid_getBlkBoundBox
  use Grid_data, ONLY : gr_globalDomain, gr_meshMe
  
  implicit none

#include "Flash.h"
#include "constants.h"
  
  integer, intent(in)           :: fileUnit ! Output to file
  logical, intent(inout)        :: perfect  ! Flag to indicate errors
  integer :: blkCount
  integer, dimension(MAXBLOCKS) :: blkList

  real, pointer, dimension(:,:,:,:) :: ctrBlk, scrBlk
  integer :: blockID, b,i,s,e
  real :: err, tolerance, blkErr, f
  real, dimension(1:MDIM) :: pos, left, right, diff
  integer, dimension(MDIM) :: ind
  integer :: procID, blkID
  real,dimension(LOW:HIGH,MDIM) :: bndbox
  logical :: tempPerfect

#ifdef FLASH_GRID_PARAMESH
  !The error is much higher when guard cell filling on an AMR mesh.
  tolerance = 0.01
#else
  tolerance = 0.000000000001
#endif

  call Grid_fillGuardCells(CENTER,ALLDIR)

  call Grid_getListOfBlocks(LEAF,blkList,blkCount)

  err = 0.0

  do b = 1,blkCount
     blockID=blkList(b)
     call Grid_getBlkPtr(blockID, ctrBlk, CENTER)
     call Grid_getBlkPtr(blockID, scrBlk, SCRATCH_CTR)
     blkErr = maxval(abs(scrBlk(1,:,:,:)-ctrBlk(1,:,:,:)))
     err = max(err,blkErr)
     call Grid_releaseBlkPtr(blockID,ctrBlk,CENTER)
     call Grid_releaseBlkPtr(blockID,scrBlk,SCRATCH_CTR)
  end do
  
  if(err>tolerance) then
     perfect=.false.
     write(fileUnit,*)'the maximum error in guardcell fill is ',err, 'the tolerance is', tolerance
  end if


  perfect=.true.
  s=1;e=7
  do b = s,e
     f=b
     pos(1:NDIM)=gr_globalDomain(LOW,1:NDIM)+&
          (gr_globalDomain(HIGH,1:NDIM)-gr_globalDomain(LOW,1:NDIM))/f
     blkID=0
     procID=0
     call Grid_getBlkIDFromPos(pos,blkList,blkCount,blkID, procID)
     if(blkID==NONEXISTENT)then
        perfect=.false.
        write(fileUnit,*)'error in finding pos, no block found'
     else if(blkID < 1 .OR. procID < 0) then
        perfect=.false.
        write(fileUnit,*)'error in finding pos, invalid (blkID,procID)=',blkID,procID
     else
        if(procID==gr_meshMe) then
           if(blkID>0) then
              call Grid_getBlkBoundBox(blkID,bndbox)
              do i = 1,NDIM
                 perfect = perfect.and.(bndbox(HIGH,i).ge.pos(i))
                 perfect = perfect.and.(bndbox(LOW,i).le.pos(i))
              end do
              write(fileUnit,1)b, procID,blkID, pos(1:NDIM),bndbox(LOW,1:NDIM),bndbox(HIGH,1:NDIM)
1             format(3I5,9F10.4)
           end if
        end if
     end if
  end do

end subroutine Grid_unitTest
