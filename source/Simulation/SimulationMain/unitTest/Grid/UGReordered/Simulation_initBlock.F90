!!****if* source/Simulation/SimulationMain/unitTest/Grid/UGReordered/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(blockId)
!!
!!
!!
!! DESCRIPTION
!!  Apply initial condition for Reordered UG unit test
!!
!! 
!! ARGUMENTS
!!
!!  blockId          the number of the block to update
!!
!!
!! PARAMETERS
!! 
!!
!!***

subroutine Simulation_initBlock(blockId)
  use Grid_interface, ONLY : Grid_getGlobalIndexLimits, &
    Grid_getBlkIndexLimits, Grid_getBlkCornerID, Grid_putPointData

!!  use Simulation_data

#include "constants.h"
#include "Flash.h"

  
  implicit none

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)
  
  integer, intent(in) :: blockId
  

  integer :: i, j, k, i1, var
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM) :: axis, globalRange
  integer, dimension(MDIM) :: startIndex,stride
  integer :: xb,xe,yb,ye,zb,ze

  real :: tempZ,tempY,tempX,temp
  real :: pi, xPi, yPi, zPi


  call Grid_getGlobalIndexLimits(globalRange)

  pi = 8.0*atan(1.0)
  xPi = pi/globalRange(IAXIS)
  yPi = pi/globalRange(JAXIS)
  zPi = pi/globalRange(KAXIS)

  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  call Grid_getBlkCornerID(blockId,startIndex,stride)
!!$  xb = 1
!!$  xe = xb + blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)
!!$
!!$  yb = 1
!!$  ye = yb + blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS) 
!!$
!!$  zb = 1
!!$  ze = zb + blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)
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
!           temp = (var-1) * 4096 + (k-1)*256 + (j-1) * 16 + (i-1)           

!           call Grid_putPointData(blockID, CENTER, var, EXTERIOR, axis, temp)
           call Grid_putPointData(blockID, CENTER, var, INTERIOR, axis, temp)
        end do

        axis(IAXIS) = axis(IAXIS)+1
      enddo
      axis(JAXIS) = axis(JAXIS)+1
    enddo
    axis(KAXIS) = axis(KAXIS) +1
  enddo
  
  return
end subroutine Simulation_initBlock
