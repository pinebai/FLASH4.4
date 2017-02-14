!!****if* source/physics/Diffuse/DiffuseMain/CG/diff_computeAX
!!
!!  NAME 
!!
!!  diff_computeAX
!!
!!  SYNOPSIS
!!
!!  call diff_computeAX (integer(IN)           :: blkLimits,
!!                  integer(IN)           :: blkLimitsGC,
!!                  real(IN)              :: X    ,
!!                  real(OUT)             :: AX ,
!!                  real(IN)              :: iFactorA    ,
!!                  real(IN)   , OPTIONAL :: iFactorC ,
!!                  real(IN)   , OPTIONAL :: iFactorD ,
!!                  real(IN)              :: dt,
!!                  real(IN)              :: theta)
!!
!!
!!  DESCRIPTION 
!!      This routine computes AX, Jacobian free, as we do not have to store AX.
!!      The elements of A are built from FD of general diffusion/conduction equation      
!!
!!      A * dV/dt = d/dx(B*dV/dx) + d/dy(B*dV/dy) + d/dx(B*dV/dz) + C*V + D
!!
!!
!! ARGUMENTS
!!      
!!      
!!
!! SIDE EFFECTS
!!      
!!  
!! NOTES:
!!  
!!
!!***

!!REORDER(4): solnData
subroutine diff_computeAX (blockID, blkLimits, blkLimitsGC, iVar, iFactorA, dt, theta, AX,iFactorC, iFactorD)
  
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &                                
       Grid_getDeltas, Grid_getBlkData
  use Timers_interface, ONLY : Timers_start, Timers_stop


  
  
  implicit none 
  
#include "Flash.h"
#include "constants.h"
  
  !!------------------------------------------------------------------------------------------
  integer, intent(IN) :: blockID
  integer, dimension(2,MDIM),intent(IN) :: blkLimits
  integer, dimension(2,MDIM),intent(IN) :: blkLimitsGC
  integer,intent(IN) :: iVar
  integer,intent(IN) :: iFactorA
  real,intent(IN):: dt
  real,intent(IN):: theta   
  real,intent(OUT) :: AX (blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))
  integer,intent(IN),OPTIONAL :: iFactorC
  integer,intent(IN),OPTIONAL :: iFactorD

  !!-------------------------------------------------------------------------------------------
  integer               :: i,j,k  
  real, dimension(MDIM) :: del
  real,POINTER,DIMENSION(:,:,:,:) :: scrch_Ctr
  real,POINTER,DIMENSION(:,:,:,:) :: solnData
  real, allocatable :: cellVolumes(:,:,:)
  
  call Timers_start("diff_computeAX")    
    
  call Grid_getBlkPtr(blockID,solnData)
  call Grid_getBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)            
  call Grid_getDeltas(blockID,del)  

  allocate(cellVolumes(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))
  
  
  call Grid_getBlkData(blockID, CELL_VOLUME, 0, EXTERIOR, &
       blkLimitsGC(LOW,:), cellVolumes(:,:,:), blkLimitsGC(HIGH,:))
  
  
  AX = 0.0
  
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)                      
           AX(i,j,k) = theta*dt*(scrch_Ctr(DFLX_SCRATCH_CENTER_VAR,i+1,j,k)-scrch_Ctr(DFLX_SCRATCH_CENTER_VAR,i,j,k))
#if NDIM >= 2
           AX(i,j,k) = AX(i,j,k) + theta*dt*(scrch_Ctr(DFLY_SCRATCH_CENTER_VAR,i,j+1,k)-scrch_Ctr(DFLY_SCRATCH_CENTER_VAR,i,j,k))           
#if NDIM == 3           
           AX(i,j,k) = AX(i,j,k) + theta*dt*(scrch_Ctr(DFLZ_SCRATCH_CENTER_VAR,i,j,k+1)-scrch_Ctr(DFLZ_SCRATCH_CENTER_VAR,i,j,k))
#endif           
#endif                
           AX(i,j,k) = AX(i,j,k) + solnData(iFactorA,i,j,k)*cellVolumes(i,j,k)* solnData(iVar,i,j,k)
           
        enddo
     end do
  end do
  
  !!-----------------------------------------------------------------------
  !!    ADD FACTORS IF PRESENT.
  !!-----------------------------------------------------------------------  
  
!!$  if (present(iFactorC)) then     
!!$     AX(:,:,:) = AX(:,:,:) + dt*solnData(iVar,:,:,:)*solnData(iFactorC,:,:,:)*cellVolumes(:,:,:)     
!!$  end if
!!$  
!!$  if (present(iFactorD)) then     
!!$     AX(:,:,:) = AX(:,:,:) + dt*solnData(iFactorD,:,:,:)*cellVolumes(:,:,:)     
!!$  end if
  
  call Grid_releaseBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)
  
  call Grid_releaseBlkPtr(blockID,solndata)
  
  deallocate(cellVolumes)
  
  call Timers_stop("diff_computeAX")       
  
end subroutine diff_computeAX
