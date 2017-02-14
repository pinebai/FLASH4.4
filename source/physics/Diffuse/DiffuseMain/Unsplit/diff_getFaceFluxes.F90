!!****if* source/physics/Diffuse/DiffuseMain/Unsplit/diff_getFaceFluxes
!!
!!  NAME 
!!
!!  diff_computeAX
!!
!!  SYNOPSIS
!!
!!  diff_computeAX (integer, intent(IN)           :: blkLimits
!!                  integer, intent(IN)           :: blkLimitsGC
!!                  real,    intent(IN)           :: X    
!!                  real,    intent(OUT)          :: AX 
!!                  real,    intent(IN)           :: iFactorB  
!!                  real,    intent(IN)           :: iFactorA    
!!                  real,    intent(IN), OPTIONAL :: iFactorC 
!!                  real,    intent(IN)           :: del
!!                  real,    intent(IN)           :: dt
!!                  real,    intent(IN)           :: theta)
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
!!  blkLimits - integer index limits for the interior of the block
!!  blkLimitsGC -integer index limits for the whole block
!!  X - input vector
!!  AX - computed value to return
!!  iFactorA -    
!!  iFactorB -    
!!  iFactorC -    
!!  del -
!!  dt - 
!!  theta - 
!!      
!!
!! SIDE EFFECTS
!!      
!!  
!! NOTES:
!!  
!!
!!***

!!REORDER(4): solnVec
subroutine diff_getFaceFluxes (blockCount, blockList, iVar, iFactorB, mode, scalefactor) 
  
  
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
  integer,intent(IN) :: mode
  real,   intent(IN) :: scalefactor
  
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
  
  if (scalefactor .NE. 0.0) then

     do lb = 1, blockCount     
     
        blockID = blockList(lb)      
     
        call Grid_getBlkPtr(blockID,solnData)   

        call Grid_getBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)          

        call Grid_getDeltas(blockID,del)     

        call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)     

        dx = del(IAXIS)
        dy = del(JAXIS)
        dz = del(KAXIS)

        !!-----------------------------------------------------------------------
        !!     1. COMPUTE THERMAL FLUXES.
        !!-----------------------------------------------------------------------

        allocate(areaLeft(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
             blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
             blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))

        call Grid_getBlkData(blockID, CELL_FACEAREA, ILO_FACE, EXTERIOR, &
             blkLimitsGC(LOW,:), areaLeft(:,:,:), blkLimitsGC(HIGH,:))

        !! Zero the fluxes.
        if (mode == 0) then
           scrch_Ctr(DFLX_SCRATCH_CENTER_VAR,:,:,:) = 0.0
#if NDIM >= 2
           scrch_Ctr(DFLY_SCRATCH_CENTER_VAR,:,:,:) = 0.0
#if NDIM == 3
           scrch_Ctr(DFLZ_SCRATCH_CENTER_VAR,:,:,:) = 0.0
#endif
#endif
        end if

        do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)+1

                 avgFacB = 0.5*(solnData(iFactorB,i,j,k)+solnData(iFactorB,i-1,j,k))*scalefactor

                 scrch_Ctr(DFLX_SCRATCH_CENTER_VAR,i,j,k) = scrch_Ctr(DFLX_SCRATCH_CENTER_VAR,i,j,k) - &
                      avgFacB*(solnData(iVar,i,j,k)-solnData(iVar,i-1,j,k))*areaLeft(i,j,k)/dx
              enddo
           end do
        end do


#if NDIM >= 2
        call Grid_getBlkData(blockID, CELL_FACEAREA, JLO_FACE, EXTERIOR, &
             blkLimitsGC(LOW,:), areaLeft(:,:,:), blkLimitsGC(HIGH,:))

        do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)+1
              do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

                 avgFacB = 0.5*(solnData(iFactorB,i,j,k)+solnData(iFactorB,i,j-1,k))*scalefactor              

                 scrch_Ctr(DFLY_SCRATCH_CENTER_VAR,i,j,k) = scrch_Ctr(DFLY_SCRATCH_CENTER_VAR,i,j,k) - &
                      avgFacB*(solnData(iVar,i,j,k)-solnData(iVar,i,j-1,k))*areaLeft(i,j,k)/dy
              enddo
           end do
        end do


#if NDIM == 3

        call Grid_getBlkData(blockID, CELL_FACEAREA, KLO_FACE, EXTERIOR, &
             blkLimitsGC(LOW,:), areaLeft(:,:,:), blkLimitsGC(HIGH,:))

        do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)+1
           do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

                 avgFacB = 0.5*(solnData(iFactorB,i,j,k)+solnData(iFactorB,i,j,k-1))*scalefactor              

                 scrch_Ctr(DFLZ_SCRATCH_CENTER_VAR,i,j,k) = scrch_Ctr(DFLZ_SCRATCH_CENTER_VAR,i,j,k) - & 
                      avgFacB*(solnData(iVar,i,j,k)-solnData(iVar,i,j,k-1))*areaLeft(i,j,k)/dz              
              enddo
           end do
        end do

#endif     
#endif    

        deallocate (areaLeft)     

        call Grid_releaseBlkPtr(blockID,solnData)   

        call Grid_releaseBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)          

     end do

  else if (mode==0) then                           ! and scalefactor==0

     ! much simplified loop
     do lb = 1, blockCount     
        blockID = blockList(lb)      
        call Grid_getBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)          
        !!-----------------------------------------------------------------------
        !!     ZERO OUT SCRATCH VARS FOR THERMAL FLUXES.
        !!-----------------------------------------------------------------------
        scrch_Ctr(DFLX_SCRATCH_CENTER_VAR,:,:,:) = 0.0
#if NDIM >= 2
        scrch_Ctr(DFLY_SCRATCH_CENTER_VAR,:,:,:) = 0.0
#if NDIM == 3
        scrch_Ctr(DFLZ_SCRATCH_CENTER_VAR,:,:,:) = 0.0
#endif
#endif
        call Grid_releaseBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)          
     end do

  end if

  !!-----------------------------------------------------------------------
  !!     2. DO FLUX CONSERVATION IF NEEDED
  !!-----------------------------------------------------------------------
  
#ifdef FLASH_GRID_PARAMESH
  
  if (mode == 1) then
     

     do lb = 1, blockCount     
        
        blockID = blockList(lb)      
        
        call Grid_getBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)          
        
        call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
        
        allocate(xflux(NFLUXES,blkLimitsGC(HIGH,IAXIS),   &
             blkLimitsGC(HIGH,JAXIS),blkLimitsGC(HIGH,KAXIS)))
     
        allocate(yflux(NFLUXES,blkLimitsGC(HIGH,IAXIS),   &
             blkLimitsGC(HIGH,JAXIS),blkLimitsGC(HIGH,KAXIS)))
        
        allocate(zflux(NFLUXES,blkLimitsGC(HIGH,IAXIS),   &
             blkLimitsGC(HIGH,JAXIS),blkLimitsGC(HIGH,KAXIS)))        
        
        xflux(1,:,:,:) =  scrch_Ctr(DFLX_SCRATCH_CENTER_VAR,:,:,:)
        call Grid_putFluxData(blockID,IAXIS,xflux,blkLimitsGC(HIGH,:))  
        
#if NDIM >= 2
        yflux(1,:,:,:) =  scrch_Ctr(DFLY_SCRATCH_CENTER_VAR,:,:,:)
        call Grid_putFluxData(blockID,JAXIS,yflux,blkLimitsGC(HIGH,:))
        
#if NDIM ==3
        zflux(1,:,:,:) =  scrch_Ctr(DFLZ_SCRATCH_CENTER_VAR,:,:,:)             
        call Grid_putFluxData(blockID,KAXIS,zflux,blkLimitsGC(HIGH,:)) 
#endif
        
#endif     
        deallocate (xflux)
        deallocate (yflux)
        deallocate (zflux)
        
        call Grid_releaseBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)     
        
     end do
     
     call Grid_conserveFluxes(ALLDIR, 0)  
     
     
     do lb = 1, blockCount     
        
        blockID = blockList(lb)      
        
        call Grid_getBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)          
        
        call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
     
        allocate(xflux(NFLUXES,blkLimitsGC(HIGH,IAXIS), &
             blkLimitsGC(HIGH,JAXIS),blkLimitsGC(HIGH,KAXIS)))
        allocate(yflux(NFLUXES,blkLimitsGC(HIGH,IAXIS), &
             blkLimitsGC(HIGH,JAXIS),blkLimitsGC(HIGH,KAXIS)))
        allocate(zflux(NFLUXES,blkLimitsGC(HIGH,IAXIS), &
             blkLimitsGC(HIGH,JAXIS),blkLimitsGC(HIGH,KAXIS)))
        
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
     
  end if
     
#endif  
  
  
  call Timers_stop("diff_computeFluxes")       
  
end subroutine diff_getFaceFluxes
