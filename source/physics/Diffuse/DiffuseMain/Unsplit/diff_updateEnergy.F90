!!****if* source/physics/Diffuse/DiffuseMain/Unsplit/diff_updateEnergy
!!
!!  NAME 
!!
!!  diff_updateEnergy
!!
!!  SYNOPSIS
!!
!!  call diff_updateEnergy (integer, intent(IN)   :: blockCount,
!!                          integer, intent(IN)   :: blockList(blockCount),
!!                          real,    intent(IN)   :: dt)
!!
!!
!!  DESCRIPTION 
!!
!!
!! ARGUMENTS
!!  dt - 
!!
!! SIDE EFFECTS
!!      
!!  
!!  
!!
!!***

!!REORDER(4): solnVec
subroutine diff_updateEnergy (blockCount, blockList, dt) 
  
  
  use Grid_interface,   ONLY : Grid_getBlkPtr,     &
                               Grid_releaseBlkPtr, &                                
                               Grid_getDeltas,     &
                               Grid_getBlkData,    &
                               Grid_getBlkIndexLimits, &
                               Grid_putFluxData, &
                               Grid_getFluxData
  use Timers_interface, ONLY : Timers_start, &
                               Timers_stop
  
  implicit none 
  
#include "Flash.h"
#include "constants.h"

  integer,intent(IN) :: blockCount
  integer,intent(IN) :: blockList(blockCount)
  real, intent(IN)   :: dt
  integer :: i,j,k  
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  integer :: blkLimits  (2,MDIM)
  integer :: blkLimitsGC(2,MDIM)
  integer :: blockID, lb   
  real, pointer, dimension(:,:,:,:) :: scrch_Ctr
  real, allocatable :: cellVolumes(:,:,:)
  real :: change
  integer :: datasize(MDIM)

  
  do lb = 1, blockCount
     
     blockID = blockList(lb)
     
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     
     allocate(cellVolumes(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)))
     
     datasize(1:MDIM)= blkLimits(HIGH,1:MDIM)-blkLimits(LOW,1:MDIM)+1        
     
     call Grid_getBlkData(blockID, CELL_VOLUME, 0, EXTERIOR,          &
          blkLimits(LOW,:), cellVolumes,datasize)
     
     call Grid_getBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)
     
     call Grid_getBlkPtr(blockID, solnVec)          
     
     
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              
              change = (scrch_Ctr(DFLX_SCRATCH_CENTER_VAR,i+1,j,k) - &
                   scrch_Ctr(DFLX_SCRATCH_CENTER_VAR,i,j,k))
              
#if NDIM >= 2
              change = change + (scrch_Ctr(DFLY_SCRATCH_CENTER_VAR,i,j+1,k) - & 
                   scrch_Ctr(DFLY_SCRATCH_CENTER_VAR,i,j,k))
#if NDIM == 3           
              change = change + (scrch_Ctr(DFLZ_SCRATCH_CENTER_VAR,i,j,k+1) - & 
                   scrch_Ctr(DFLZ_SCRATCH_CENTER_VAR,i,j,k))
#endif        
              
#endif               
              
#ifndef TELE_VAR
              solnVec(EINT_VAR,i,j,k) = (solnVec(DENS_VAR,i,j,k)*solnVec(EINT_VAR,i,j,k) - & 
                   (change/cellVolumes(i,j,k)))/solnVec(DENS_VAR,i,j,k)
              solnVec(ENER_VAR,i,j,k) = (solnVec(DENS_VAR,i,j,k)*solnVec(ENER_VAR,i,j,k) - & 
                   (change/cellVolumes(i,j,k)))/solnVec(DENS_VAR,i,j,k)
#else              
              solnVec(EELE_VAR,i,j,k) = (solnVec(DENS_VAR,i,j,k)*solnVec(EELE_VAR,i,j,k) - & 
                   (dt*change/cellVolumes(i,j,k)))/solnVec(DENS_VAR,i,j,k)                           
#endif
              
           end do
        end do
     end do
     
     deallocate (cellVolumes)
     
     call Grid_releaseBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)
     
     call Grid_releaseBlkPtr(blockID,solnVec)
     
  end do
  
end subroutine diff_updateEnergy
