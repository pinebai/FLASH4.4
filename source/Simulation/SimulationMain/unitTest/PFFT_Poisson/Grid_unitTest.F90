!!****if* source/Simulation/SimulationMain/unitTest/PFFT_Poisson/Grid_unitTest
!!
!! NAME
!!
!!  Grid_unitTest
!!
!! SYNOPSIS
!!
!!  Grid_unitTest(integer, intent(in):: fileUnit,
!!                logical, intent(inout)::perfect  )
!!
!! DESCRIPTION
!!
!!  This unit test exercises the data accessing functions of the Grid unit
!!  The routine has direct access to all the mesh data structures such as 
!!  "unk", "facevarx" etc. It uses the Grid_getBlk/Point/RowData functions 
!!  to fetch some or all of the block data, and then compares it with
!!  the corresponding section of the appropriate array.
!!
!! ARGUMENTS
!!
!!  fileUnit - open f90 write unit
!!  perfect - returns a true if the test passed, false otherwise
!!
!! NOTES
!!
!!***

subroutine Grid_unitTest(fileUnit,perfect)

  use Grid_interface, ONLY : Grid_solvePoisson,Grid_getBlkIndexLimits,&
       Grid_getBlkPtr,Grid_releaseBlkPtr, Grid_getListOfBlocks

  implicit none

#include "Flash.h"
#include "constants.h"
  
  integer, intent(in)           :: fileUnit ! Output to file
  logical, intent(inout)        :: perfect  ! Flag to indicate errors
  
  real, pointer, dimension(:,:,:,:) :: solnData
  
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  integer, dimension(2*MDIM) :: bcTypes
  real,dimension(2,2*MDIM) :: bcValues
  real :: err, cumulativeErr,poisfact
  integer,dimension(MAXBLOCKS) :: blkList
  integer :: blockID,blkCount,lb,i,j,k

  cumulativeErr = 0.0
  poisfact = 4.0 * PI * PI
  call Grid_solvePoisson(PFFT_VAR,DENS_VAR,bcTypes,bcValues,poisfact)
  
  call Grid_getListOfBlocks(LEAF,blkList,blkCount)


  if (NDIM .EQ. 1) then
    poisfact = -4.0
  else
     poisfact=-13.0
  endif
  do lb = 1,blkCount
     blockID=blkList(lb)

     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blockID,solnData,CENTER)

     err=maxval(abs(&
          solnData(DENS_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
                            blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
                            blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))-&
          solnData(PFFT_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
                            blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
                            blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))*&
                            poisfact))
     cumulativeErr=max(cumulativeErr,err)

     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              solnData(DIFF_VAR,i,j,k)=solnData(DENS_VAR,i,j,k)-poisfact*solnData(PFFT_VAR,i,j,k)
           end do
        end do
     end do
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
  end do
  
  perfect=(abs(cumulativeErr) < 1e-05)
  
  print*,"the result is ", perfect,cumulativeErr,err
  
  return
 
end subroutine Grid_unitTest
