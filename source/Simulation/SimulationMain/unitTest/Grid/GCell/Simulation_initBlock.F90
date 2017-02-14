!!****if* source/Simulation/SimulationMain/unitTest/Grid/GCell/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer (IN) ::blockId)
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes the Grid with a sinusoidal function, volume average of cos(x)*cos(y)*cos(z)
!! 
!! ARGUMENTS
!!
!!  blockId -          the blockId to update
!!
!!
!!***

subroutine Simulation_initBlock(blockId)
  use Grid_interface, ONLY : Grid_getGlobalIndexLimits, &
    Grid_getBlkIndexLimits, Grid_getBlkCornerID, Grid_putPointData, &
    Grid_getCellCoords, Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getDeltas

#include "constants.h"
#include "Flash.h"

  
  implicit none

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)
  
  integer, intent(in) :: blockId

  real :: twopi, dx, dy, dz

  integer :: i, j, k
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC

  integer, dimension(MDIM) :: blkSize
  real,pointer,dimension(:,:,:,:)::solnData
  real,pointer,dimension(:,:,:,:):: gridVarData

  real,dimension(:),allocatable :: xleft, xright, yleft, yright, zright, zleft

  real, dimension(MDIM) :: delta
  
  twopi = PI*2.0

  
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  blkSize=blkLimitsGC(HIGH,:)-blkLimitsGC(LOW,:)+1

 
  allocate(xleft(blkSize(IAXIS)))
  allocate(xright(blkSize(IAXIS)))
  allocate(yleft(blkSize(JAXIS)))
  allocate(yright(blkSize(JAXIS)))
  allocate(zleft(blkSize(KAXIS)))
  allocate(zright(blkSize(KAXIS)))
  call Grid_getCellCoords(IAXIS,blockID,LEFT_EDGE,.true.,xleft,blkSize(IAXIS))
  call Grid_getCellCoords(IAXIS,blockID,RIGHT_EDGE,.true.,xright,blkSize(IAXIS))
  call Grid_getCellCoords(JAXIS,blockID,LEFT_EDGE,.true.,yleft,blkSize(JAXIS))
  call Grid_getCellCoords(JAXIS,blockID,RIGHT_EDGE,.true.,yright,blkSize(JAXIS))
  call Grid_getCellCoords(KAXIS,blockID,LEFT_EDGE,.true.,zleft,blkSize(KAXIS))
  call Grid_getCellCoords(KAXIS,blockID,RIGHT_EDGE,.true.,zright,blkSize(KAXIS))

 call Grid_getDeltas(blockID, delta)
 dx = delta(IAXIS)
 dy = delta(JAXIS)
 dz = delta(KAXIS)


  call Grid_getBlkPtr(blockID,solnData,CENTER)
  do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           solnData(:,i,j,k)=(1/(dx*twopi))*(sin(twopi*xright(i)) - sin(twopi*xleft(i)))
           if(NDIM>1)&
                solnData(:,i,j,k)=solnData(:,i,j,k)*&
                (1/(dy*twopi))*(sin(twopi*yright(j)) - sin(twopi*yleft(j)))
           if(NDIM>2)&
                solnData(:,i,j,k)=solnData(:,i,j,k)*&
                (1/(dz*twopi))*(sin(twopi*zright(k)) - sin(twopi*zleft(k)))
        end do
     end do
  end do
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)
  call Grid_getBlkPtr(blockID,gridVarData,SCRATCH_CTR)
  do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
           gridVarData(:,i,j,k)=(1/(dx*twopi))*(sin(twopi*xright(i)) - sin(twopi*xleft(i)))
           if(NDIM>1)&
                gridVarData(:,i,j,k)=gridVarData(:,i,j,k)*&
                (1/(dy*twopi))*(sin(twopi*yright(j)) - sin(twopi*yleft(j)))
           if(NDIM>2)&
                gridVarData(:,i,j,k)=gridVarData(:,i,j,k)*&
                (1/(dz*twopi))*(sin(twopi*zright(k)) - sin(twopi*zleft(k)))
        end do
     end do
  end do
  call Grid_releaseBlkPtr(blockID,gridVarData,SCRATCH_CTR)
  deallocate(xright)
  deallocate(yright)
  deallocate(zright)
  deallocate(xleft)
  deallocate(yleft)
  deallocate(zleft)

  return
end subroutine Simulation_initBlock
