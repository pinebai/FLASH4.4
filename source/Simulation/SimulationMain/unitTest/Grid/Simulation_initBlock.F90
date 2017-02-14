!!****if* source/Simulation/SimulationMain/unitTest/Grid/Simulation_initBlock
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
!!  Initializes the Grid with a composit number which is a combination
!!  of the block number and the indices of the cell
!! 
!! ARGUMENTS
!!
!!  blockId -          the blockId to update
!!  
!!
!!
!!***

subroutine Simulation_initBlock(blockId)
  use Grid_interface, ONLY :  Grid_putPointData, Grid_getBlkIndexLimits
  use Grid_data, ONLY : gr_meshME
#include "constants.h"
#include "Flash.h"

  
  implicit none

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)
  
  integer, intent(in) :: blockId
  


  integer :: i, j, k, var
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM) :: grd,axis
  
  real :: blk,kind,jind,iind,val

  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  grd(:)=blkLimits(LOW,:)-blkLimitsGC(LOW,:)

  blk=(gr_MeshMe*100.0+blockID)*100000000.0
  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     kind = k*1000000.0
     axis(KAXIS)=k
    do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
       jind = j*10000.0
       axis(JAXIS)=j
       do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
          iind = i*100.0
          axis(IAXIS)=i
          do var=UNK_VARS_BEGIN,UNK_VARS_END
             val = blk+kind+jind+iind+var
             call Grid_putPointData(blockID, CENTER, var, EXTERIOR, axis, val)
          end do
       enddo
    enddo
 enddo
  
 return
end subroutine Simulation_initBlock
