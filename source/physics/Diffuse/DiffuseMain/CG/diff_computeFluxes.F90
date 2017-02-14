!!****if* source/physics/Diffuse/DiffuseMain/CG/diff_computeFluxes
!!
!!  NAME 
!!
!!  diff_computeFluxes
!!
!!  SYNOPSIS
!!
!!  call diff_computeFluxes (integer(IN)           :: blockCount,
!!                              real(IN)           :: iFactorB)
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
subroutine diff_computeFluxes (blockCount,blockList,iVar,iFactorB) 
  
  
  use Grid_interface,   ONLY : Grid_getBlkPtr,     &
                               Grid_releaseBlkPtr, &                                
                               Grid_getDeltas,     &
                               Grid_getBlkData,    &
                               Grid_getBlkIndexLimits, &
                               Grid_putFluxData, &
                               Grid_getFluxData
  use Timers_interface, ONLY : Timers_start, &
                               Timers_stop

  use Diffuse_data, ONLY     : diff_geometricMeanDiff
  
  implicit none 
  
#include "Flash.h"
#include "constants.h"

  integer,intent(IN) :: blockCount
  integer,intent(IN) :: blockList(blockCount)
  integer,intent(IN) :: iVar 
  integer,intent(IN) :: iFactorB
  
  integer :: i,j,k  
  real, POINTER, DIMENSION(:,:,:,:) :: solnData
  real, dimension(MDIM) :: del
  integer :: datasizeGC(MDIM)
  real, allocatable :: areaLeft(:,:,:)
  real :: avgFacB
  integer :: blkLimits  (2,MDIM)
  integer :: blkLimitsGC(2,MDIM)
  real :: dx,dy,dz
  integer :: blockID, lb
  
  real, allocatable :: xflux(:,:,:,:)
  real, allocatable :: yflux(:,:,:,:)
  real, allocatable :: zflux(:,:,:,:)
  
  real, pointer, dimension(:,:,:,:) :: scrch_Ctr
  
  call Timers_start("diff_computeFluxes")
  
  do lb = 1, blockCount     
     blockID = blockList(lb)      
     call Grid_getBlkPtr(blockID,solnData)   
     call Grid_getBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)          
     call Grid_getDeltas(blockID,del)     
     call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
     
     allocate(xflux(NFLUXES,blkLimitsGC(HIGH,IAXIS),blkLimitsGC(HIGH,JAXIS),blkLimitsGC(HIGH,KAXIS)))
     allocate(yflux(NFLUXES,blkLimitsGC(HIGH,IAXIS),blkLimitsGC(HIGH,JAXIS),blkLimitsGC(HIGH,KAXIS)))
     allocate(zflux(NFLUXES,blkLimitsGC(HIGH,IAXIS),blkLimitsGC(HIGH,JAXIS),blkLimitsGC(HIGH,KAXIS)))
  
     dx = del(IAXIS)
     dy = del(JAXIS)
     dz = del(KAXIS)

     xflux = 0.0
     yflux = 0.0
     zflux = 0.0
  
     !!-----------------------------------------------------------------------
     !!     1. COMPUTE THERMAL FLUXES.
     !!-----------------------------------------------------------------------
     
     allocate(areaLeft(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))

     call Grid_getBlkData(blockID, CELL_FACEAREA, ILO_FACE, EXTERIOR, &
          blkLimitsGC(LOW,:), areaLeft(:,:,:), blkLimitsGC(HIGH,:))

     
     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)+1
              
              if (diff_geometricMeanDiff) then
                 avgFacB = 2.0*solnData(iFactorB,i,j,k)*solnData(iFactorB,i-1,j,k)/ &
                      (solnData(iFactorB,i,j,k)+solnData(iFactorB,i-1,j,k))
              else              
                 avgFacB = 0.5*(solnData(iFactorB,i,j,k)+solnData(iFactorB,i-1,j,k))              
              endif
              scrch_Ctr(DFLX_SCRATCH_CENTER_VAR,i,j,k) = -avgFacB*(solnData(iVar,i,j,k)-solnData(iVar,i-1,j,k))*areaLeft(i,j,k)/dx
           enddo
        end do
     end do     
     
     xflux(1,:,:,:) =  scrch_Ctr(DFLX_SCRATCH_CENTER_VAR,:,:,:)
     
     call Grid_putFluxData(blockID,IAXIS,xflux,blkLimitsGC(HIGH,:))          
     
     
#if NDIM >= 2
     call Grid_getBlkData(blockID, CELL_FACEAREA, JLO_FACE, EXTERIOR, &
          blkLimitsGC(LOW,:), areaLeft(:,:,:), blkLimitsGC(HIGH,:))
     
     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)+1
           do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              
              if (diff_geometricMeanDiff) then
                 avgFacB = 2.0*solnData(iFactorB,i,j,k)*solnData(iFactorB,i,j-1,k)/ &
                      (solnData(iFactorB,i,j,k)+solnData(iFactorB,i,j-1,k))
              else              
                 avgFacB = 0.5*(solnData(iFactorB,i,j,k)+solnData(iFactorB,i,j-1,k))              
              endif
              scrch_Ctr(DFLY_SCRATCH_CENTER_VAR,i,j,k) = -avgFacB*(solnData(iVar,i,j,k)-solnData(iVar,i,j-1,k))*areaLeft(i,j,k)/dy
           enddo
        end do
     end do

     yflux(1,:,:,:) =  scrch_Ctr(DFLY_SCRATCH_CENTER_VAR,:,:,:)
     
     call Grid_putFluxData(blockID,JAXIS,yflux,blkLimitsGC(HIGH,:))     

#if NDIM == 3

     call Grid_getBlkData(blockID, CELL_FACEAREA, KLO_FACE, EXTERIOR, &
          blkLimitsGC(LOW,:), areaLeft(:,:,:), blkLimitsGC(HIGH,:))
     
     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)+1
        do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              
              if (diff_geometricMeanDiff) then
                 avgFacB = 2.0*solnData(iFactorB,i,j,k)*solnData(iFactorB,i,j,k-1)/ &
                      (solnData(iFactorB,i,j,k)+solnData(iFactorB,i,j,k-1))
              else              
                 avgFacB = 0.5*(solnData(iFactorB,i,j,k)+solnData(iFactorB,i,j,k-1))              
              endif
              scrch_Ctr(DFLZ_SCRATCH_CENTER_VAR,i,j,k) = -avgFacB*(solnData(iVar,i,j,k)-solnData(iVar,i,j,k-1))*areaLeft(i,j,k)/dz
           enddo
        end do
     end do
     
     zflux(1,:,:,:) =  scrch_Ctr(DFLZ_SCRATCH_CENTER_VAR,:,:,:)     
     
     call Grid_putFluxData(blockID,KAXIS,zflux,blkLimitsGC(HIGH,:))     
#endif     
#endif
     
     deallocate (xflux)
     deallocate (yflux)
     deallocate (zflux)

     deallocate (areaLeft)     
     
     call Grid_releaseBlkPtr(blockID,solnData)   
     call Grid_releaseBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)          

     
  end do
  
  !!-----------------------------------------------------------------------
  !!     2. DO FLUX CONSERVATION IF NEEDED
  !!-----------------------------------------------------------------------
  
#ifdef FLASH_GRID_PARAMESH

  call Grid_conserveFluxes(IAXIS, 0)  
#if NDIM >= 2  
  call Grid_conserveFluxes(JAXIS, 0)
#if NDIM == 3  
  call Grid_conserveFluxes(KAXIS, 0)  
#endif
#endif  
  
  do lb = 1, blockCount     
     blockID = blockList(lb)      
     call Grid_getBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)          
     call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
     
     allocate(xflux(NFLUXES,blkLimitsGC(HIGH,IAXIS),blkLimitsGC(HIGH,JAXIS),blkLimitsGC(HIGH,KAXIS)))
     allocate(yflux(NFLUXES,blkLimitsGC(HIGH,IAXIS),blkLimitsGC(HIGH,JAXIS),blkLimitsGC(HIGH,KAXIS)))
     allocate(zflux(NFLUXES,blkLimitsGC(HIGH,IAXIS),blkLimitsGC(HIGH,JAXIS),blkLimitsGC(HIGH,KAXIS)))

     xflux = 0.0
     yflux = 0.0
     zflux = 0.0
     
     call Grid_getFluxData(blockID,IAXIS,xflux,blkLimitsGC(HIGH,:))          
     
     scrch_Ctr(DFLX_SCRATCH_CENTER_VAR,blkLimits(LOW, IAXIS),:,:) = xflux(1,blkLimits(LOW,IAXIS),:,:)
     scrch_Ctr(DFLX_SCRATCH_CENTER_VAR,blkLimits(HIGH,IAXIS)+1,:,:) = xflux(1,blkLimits(HIGH,IAXIS)+1,:,:)

#if NDIM >= 2
     call Grid_getFluxData(blockID,JAXIS,yflux,blkLimitsGC(HIGH,:))          
     scrch_Ctr(DFLY_SCRATCH_CENTER_VAR,:,blkLimits(LOW, JAXIS),:) = yflux(1,:,blkLimits(LOW,JAXIS),:)
     scrch_Ctr(DFLY_SCRATCH_CENTER_VAR,:,blkLimits(HIGH,JAXIS)+1,:) = yflux(1,:,blkLimits(HIGH,JAXIS)+1,:)

#if NDIM == 3
     call Grid_getFluxData(blockID,KAXIS,zflux,blkLimitsGC(HIGH,:))          
     scrch_Ctr(DFLZ_SCRATCH_CENTER_VAR,:,:,blkLimits(LOW, KAXIS))   = zflux(1,:,:,blkLimits(LOW,KAXIS))
     scrch_Ctr(DFLZ_SCRATCH_CENTER_VAR,:,:,blkLimits(HIGH,KAXIS)+1) = zflux(1,:,:,blkLimits(HIGH,KAXIS)+1)
#endif
     
#endif               
     
     deallocate (xflux)
     deallocate (yflux)
     deallocate (zflux)
     
     call Grid_releaseBlkPtr(blockID,solnData)   
     
     call Grid_releaseBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)     
     
  end do
  
#endif

  
  
  call Grid_releaseBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)
  
  call Timers_stop("diff_computeFluxes")       
  
end subroutine
