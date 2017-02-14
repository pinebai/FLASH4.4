!!****if* source/Simulation/SimulationMain/unitTest/Grid/UG/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer (IN) ::blockId, 
!!
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes the Grid with a sinusoidal function sin(x)*cos(y)*cos(z)
!! 
!! ARGUMENTS
!!
!!  blockId -          the blockId to update
!!  
!!
!!
!!***

subroutine Simulation_initBlock(blockId)
  use Grid_interface, ONLY : Grid_getGlobalIndexLimits, &
    Grid_getBlkIndexLimits, Grid_getBlkCornerID, Grid_putPointData, &
    Grid_getCellCoords, Grid_getBlkPtr, Grid_releaseBlkPtr

#include "constants.h"
#include "Flash.h"

  
  implicit none

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)
  
  integer, intent(in) :: blockId
  

  integer, dimension(MDIM) :: axis, globalRange
  integer, dimension(MDIM) :: startIndex,stride
  integer :: xb,xe,yb,ye,zb,ze

  real :: tempZ,tempY,tempX,temp
  real :: twopi, xPi, yPi, zPi

  integer :: i, j, k, i1, var
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC

  integer, dimension(MDIM) :: blkSize
  real,pointer,dimension(:,:,:,:)::solnData

  real,dimension(:),allocatable :: xCoords, yCoords, zCoords
  
  twopi = PI*2.0


  call Grid_getGlobalIndexLimits(globalRange)

  xPi = twopi/globalRange(IAXIS)
  yPi = twopi/globalRange(JAXIS)
  zPi = twopi/globalRange(KAXIS)

  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  call Grid_getBlkCornerID(blockId,startIndex,stride)
  xb = startIndex(IAXIS)
  xe = xb + blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)

  yb = startIndex(JAXIS)
  ye = yb + blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)

  zb = startIndex(KAXIS)
  ze = zb + blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)

  axis(KAXIS) = 1
  do k = zb,ze
    ! alternate initialization values
    ! tempZ = 100*k
    tempZ = cos((k-1)*zPi)
    axis(JAXIS) = 1
    do j = yb,ye
      ! tempY = (tempZ+j)*100
      tempY = cos((j-1)*yPi)
      axis(IAXIS) = 1
      do i = xb,xe

        tempX = sin((i-1)*xPi)
        temp=tempX*tempY*tempZ
        do var=UNK_VARS_BEGIN,UNK_VARS_END

           call Grid_putPointData(blockID, CENTER, var, INTERIOR, axis, temp)
        end do

        axis(IAXIS) = axis(IAXIS)+1
      enddo
      axis(JAXIS) = axis(JAXIS)+1
    enddo
    axis(KAXIS) = axis(KAXIS) +1
  enddo
  
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,FACEX)
  blkSize=blkLimitsGC(HIGH,:)-blkLimitsGC(LOW,:)+1
  allocate(xCoords(blkSize(IAXIS)))
  allocate(yCoords(blkSize(JAXIS)))
  allocate(zCoords(blkSize(KAXIS)))
  call Grid_getBlkPtr(blockID,solnData,FACEX)
  call Grid_getCellCoords(IAXIS,blockID,FACES,.true.,xCoords,blkSize(IAXIS))
  call Grid_getCellCoords(JAXIS,blockID,CENTER,.true.,yCoords,blkSize(JAXIS))
  call Grid_getCellCoords(KAXIS,blockID,CENTER,.true.,zCoords,blkSize(KAXIS))
  do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           solnData(:,i,j,k)=cos(twopi*xCoords(i))*cos(twopi*yCoords(j))&
                *cos(twopi*zCoords(k))
        end do
     end do
  end do
  call Grid_releaseBlkPtr(blockID,solnData,FACEX)
  deallocate(xCoords)
  deallocate(yCoords)
  deallocate(zCoords)

#if (NDIM>1)
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,FACEY)
  blkSize=blkLimitsGC(HIGH,:)-blkLimitsGC(LOW,:)+1
  allocate(xCoords(blkSize(IAXIS)))
  allocate(yCoords(blkSize(JAXIS)))
  allocate(zCoords(blkSize(KAXIS)))
  call Grid_getBlkPtr(blockID,solnData,FACEY)
  call Grid_getCellCoords(IAXIS,blockID,CENTER,.true.,xCoords,blkSize(IAXIS))
  call Grid_getCellCoords(JAXIS,blockID,FACES,.true.,yCoords,blkSize(JAXIS))
  call Grid_getCellCoords(KAXIS,blockID,CENTER,.true.,zCoords,blkSize(KAXIS))
  do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           solnData(:,i,j,k)=cos(twopi*xCoords(i))*cos(twopi*yCoords(j))&
                *cos(twopi*zCoords(k))
        end do
     end do
  end do
  call Grid_releaseBlkPtr(blockID,solnData,FACEY)
  deallocate(xCoords)
  deallocate(yCoords)
  deallocate(zCoords)
#endif

#if(NDIM>2)
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,FACEZ)
  blkSize=blkLimitsGC(HIGH,:)-blkLimitsGC(LOW,:)+1
  allocate(xCoords(blkSize(IAXIS)))
  allocate(yCoords(blkSize(JAXIS)))
  allocate(zCoords(blkSize(KAXIS)))
  call Grid_getBlkPtr(blockID,solnData,FACEZ)
  call Grid_getCellCoords(IAXIS,blockID,CENTER,.true.,xCoords,blkSize(IAXIS))
  call Grid_getCellCoords(JAXIS,blockID,CENTER,.true.,yCoords,blkSize(JAXIS))
  call Grid_getCellCoords(KAXIS,blockID,FACES,.true.,zCoords,blkSize(KAXIS))
  do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           solnData(:,i,j,k)=cos(twopi*xCoords(i))*cos(twopi*yCoords(j))&
                *cos(twopi*zCoords(k))
        end do
     end do
  end do
  call Grid_releaseBlkPtr(blockID,solnData,FACEZ)
  deallocate(xCoords)
  deallocate(yCoords)
  deallocate(zCoords)
#endif
  
  return
end subroutine Simulation_initBlock
