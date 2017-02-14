!!****if* source/Simulation/SimulationMain/unitTest/PFFT_Poisson/Simulation_initBlock
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
!!  Initializes the Grid with some sort of sinusoidal function in x & y
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
       Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr,&
       Grid_getBlkCornerID

#include "constants.h"
#include "Flash.h"

  
  implicit none

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)
  
  integer, intent(in) :: blockId
  

  integer, dimension(MDIM) :: axis, globalSize,startIndex,stride

  integer :: i, j, k, ii,jj,kk
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC

  integer, dimension(MDIM) :: blkSize
  real,pointer,dimension(:,:,:,:)::solnData
  real :: xPi, yPi, zPi, twopi

  twopi = PI*2.0

  call Grid_getGlobalIndexLimits(globalSize)


  xPi = twopi/globalSize(IAXIS)
  
  yPi = twopi/globalSize(JAXIS)
  zPi = twopi/globalSize(KAXIS)

  call Grid_getBlkCornerID(blockId,startIndex,stride)

  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)

  call Grid_getBlkPtr(blockID,solnData,CENTER)

  
  do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     jj=startIndex(JAXIS)
     do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        ii=startIndex(IAXIS)
        do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           if(NDIM==1) then
              solnData(DENS_VAR,i,j,k)=-4.0*cos(2.0*xPi*(ii-1))
           else
              solnData(DENS_VAR,i,j,k)=-13.0*cos(2.0*xPi*(ii-1))*sin(3.0*yPi*(jj-1))
           end if
           ii=ii+1
        end do
        jj=jj+1
     end do
  end do

  call Grid_releaseBlkPtr(blockID,solnData,CENTER)
  return
end subroutine Simulation_initBlock
