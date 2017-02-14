!!****if* source/physics/Diffuse/DiffuseMain/CG/diff_computeblkAMat
!!
!!  NAME 
!!
!!  diff_computeblkAMat
!!
!!  SYNOPSIS
!!
!!  call diff_computeblkAMat(AA, JA, IA, UPTR, NZ, N, blockID, iFactorA, iFactorB, iFactorC)!!
!!
!!  DESCRIPTION 
!!
!!
!!
!!
!! ARGUMENTS
!!
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES
!!
!!
!!***

!!REORDER(4): solnVec

subroutine diff_computeblkAMat (blockID, NZ, N, dt, theta, AA, JA, IA, UPTR, iFactorA, iFactorB, iFactorC)
  
  use Grid_interface,    ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
                                Grid_getBlkIndexLimits, &
                                Grid_getDeltas, Grid_getBlkBC
  use Timers_interface,  ONLY : Timers_start, Timers_stop
  
  implicit none
  
#include "Flash.h"
#include "constants.h"  
  
  integer, intent (IN)  :: blockID
  integer, intent (IN)  :: NZ
  integer, intent (IN)  :: N 
  real,    intent (IN)  :: dt
  real,    intent (IN)  :: theta  
  real,    intent (OUT) :: AA  (NZ)
  integer, intent (OUT) :: JA  (NZ)
  integer, intent (OUT) :: IA  (N+1)
  integer, intent (OUT) :: UPTR(N)  
  integer, intent(IN)   :: iFactorA
  integer, intent(IN)   :: iFactorB
  integer, OPTIONAL, intent(IN) :: iFactorC  
  
  real, POINTER, DIMENSION(:,:,:,:) :: solnData
  integer, DIMENSION(2,MDIM) :: blkLimitsGC, blkLimits 
  real, dimension(MDIM) :: del   

  real    :: CDiv
  integer :: pos_ijk, ia_iter, cs_iter
    
  real   :: condiph, condimh
  real   :: condjph, condjmh
  real   :: condkph, condkmh 

  integer :: i, j, k, pos, k1, k2
  integer :: datasize(MDIM), faces(2,MDIM)
  
  integer    :: add2diag 

  real, allocatable :: cellVolumes(:,:,:)
  real, allocatable :: faceAreas  (:,:,:,:) 


    
  
  call Timers_start("diff_computeblkAMat")  

  AA = 0.0
  
  cs_iter = 1
  ia_iter = 1
  
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)    
  call Grid_getDeltas(blockID, del)
  call Grid_getBlkPtr(blockID, solnData)      

  allocate(cellVolumes(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))
  
  call Grid_getBlkData(blockID, CELL_VOLUME, 0, EXTERIOR, &
       blkLimitsGC(LOW,:), cellVolumes(:,:,:), blkLimitsGC(HIGH,:))


  
  
  allocate(faceAreas(2*NDIM,                            & 
       blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),  &
       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),  &
       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))            

  call Grid_getBlkData(blockID, CELL_FACEAREA, ILO_FACE, EXTERIOR, & 
       blkLimitsGC(LOW,:), faceAreas(ILO_FACE,:,:,:), blkLimitsGC(HIGH,:))     
  
  call Grid_getBlkData(blockID, CELL_FACEAREA, IHI_FACE, EXTERIOR, &
       blkLimitsGC(LOW,:), faceAreas(IHI_FACE,:,:,:), blkLimitsGC(HIGH,:))     
  
#if NDIM >= 2
  
  call Grid_getBlkData(blockID, CELL_FACEAREA, JLO_FACE, EXTERIOR, &
       blkLimitsGC(LOW,:), faceAreas(JLO_FACE,:,:,:), blkLimitsGC(HIGH,:))   
  
  call Grid_getBlkData(blockID, CELL_FACEAREA, JHI_FACE, EXTERIOR, &
       blkLimitsGC(LOW,:), faceAreas(JHI_FACE,:,:,:), blkLimitsGC(HIGH,:))     
  
#if NDIM == 3
  
  call Grid_getBlkData(blockID, CELL_FACEAREA, KLO_FACE, EXTERIOR, &
       blkLimitsGC(LOW,:), faceAreas(KLO_FACE,:,:,:), blkLimitsGC(HIGH,:))   
  
  call Grid_getBlkData(blockID, CELL_FACEAREA, KHI_FACE, EXTERIOR, &
       blkLimitsGC(LOW,:), faceAreas(KHI_FACE,:,:,:), blkLimitsGC(HIGH,:))  
#endif

#endif

  datasize(1:MDIM) = blkLimits(HIGH,1:MDIM)- blkLimits(LOW,1:MDIM) + 1   
  
  do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)      
     do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)                   

           add2diag = 0
           cDiv = solnData(iFactorA,i,j,k)           
           
           !! i,j,k => position in matrix            
           
           pos_ijk = (i-blkLimits(LOW,IAXIS)+1)                                +  &
                     (j-blkLimits(LOW,JAXIS))*dataSize(IAXIS)                  +  &
                     (k-blkLimits(LOW,KAXIS))*dataSize(IAXIS)*dataSize(JAXIS)


           
           condiph = 0.5*(solnData(iFactorB,i+1,j,k)+ solnData(iFactorB,i,j,k))*faceAreas(ILO_FACE,i,j,k)
           condimh = 0.5*(solnData(iFactorB,i-1,j,k)+ solnData(iFactorB,i,j,k))*faceAreas(IHI_FACE,i,j,k)
#if NDIM >= 2             
           condjph = 0.5*(solnData(iFactorB,i,j+1,k)+ solnData(iFactorB,i,j,k))*faceAreas(JLO_FACE,i,j,k)
           condjmh = 0.5*(solnData(iFactorB,i,j-1,k)+ solnData(iFactorB,i,j,k))*faceAreas(JHI_FACE,i,j,k)             
#if NDIM == 3             
           condkph = 0.5*(solnData(iFactorB,i,j,k+1)+ solnData(iFactorB,i,j,k))*faceAreas(KLO_FACE,i,j,k)
           condkmh = 0.5*(solnData(iFactorB,i,j,k-1)+ solnData(iFactorB,i,j,k))*faceAreas(KHI_FACE,i,j,k) 
#endif  
#endif
           
           IA(ia_iter) = cs_iter    
           
#if NDIM == 3
           !! i,j,k-1
           if (k /= blkLimits(LOW, KAXIS)) then      
              AA(cs_iter) = -condkmh*theta*dt/del(KAXIS)
              pos = pos_ijk - datasize(IAXIS)*datasize(JAXIS)
              JA(cs_iter) = pos
              cs_iter = cs_iter + 1  
           end if
#endif
           
#if NDIM >= 2           
           !! i,j-1,k
           if (j /= blkLimits(LOW, JAXIS)) then              
              AA(cs_iter) = -condjmh*theta*dt/del(JAXIS)
              pos = pos_ijk - datasize(IAXIS)
              JA(cs_iter) = pos
              cs_iter = cs_iter + 1


           end if
#endif
           !! i-1,j,k
           if (i /= blkLimits(LOW, IAXIS)) then              
              AA(cs_iter) = -condimh*theta*dt/del(IAXIS)
              pos = pos_ijk - 1
              JA(cs_iter) = pos              
              cs_iter = cs_iter + 1                     

           end if
           
           !! i,j,k
           AA(cs_iter) = CDiv*cellVolumes(i,j,k) 
           AA(cs_iter) =  AA(cs_iter) + (theta*(Condimh+Condiph)*dt/del(IAXIS))
#if NDIM >= 2           
           AA(cs_iter) =  AA(cs_iter) + (theta*(Condjmh+Condjph)*dt/del(JAXIS))*K2D
#if NDIM == 3
           AA(cs_iter) =  AA(cs_iter) + (theta*(Condkmh+Condkph)*dt/del(KAXIS))*K3D               
#endif
#endif
           

           
           if (present(iFactorC)) then               
              AA(cs_iter) = AA(cs_iter) + dt * solnData(iFactorC,i,j,k)*cellVolumes(i,j,k)
           end if
           
           UPTR(ia_iter) = cs_iter
           pos = pos_ijk
           JA(cs_iter) = pos            
           cs_iter = cs_iter + 1        


           
           !! i+1,j,k
           if (i /= blkLimits(HIGH, IAXIS)) then              
              AA(cs_iter) = -condiph*theta*dt/del(IAXIS)
              pos = pos_ijk + 1
              JA(cs_iter) = pos
              cs_iter = cs_iter + 1      
                            
           end if
           
#if NDIM >= 2
           !! i,j+1,k
           if (j /= blkLimits(HIGH, JAXIS)) then              
              AA(cs_iter) = -condjph*theta*dt/del(JAXIS)  
              pos = pos_ijk + datasize(IAXIS)
              JA(cs_iter) = pos
              cs_iter = cs_iter + 1   
              
           end if
#endif
           
           
#if NDIM == 3
           !! i,j,k+1
           if (k /= blkLimits(HIGH, KAXIS)) then
              AA(cs_iter) = -condkph*theta*dt/del(KAXIS)  
              pos = pos_ijk + datasize(IAXIS)*datasize(JAXIS)
              JA(cs_iter) = pos
              cs_iter = cs_iter + 1  
              
           end if
#endif                    
           
!!$           AA(IA(ia_iter):cs_iter) = AA(IA(ia_iter):cs_iter) / AA(UPTR(ia_iter))
!!$           AA(UPTR(ia_iter)) = 1.0

!!$           write (*,*) AA(IA(ia_iter):cs_iter-1)
           
           ia_iter = ia_iter + 1


              

           
        end do
     end do
  end do
  
  IA(ia_iter) = 1 + NZ
  
  call Grid_releaseBlkPtr(blockID, solnData)   
  
  
  deallocate (cellVolumes)
  deallocate (faceAreas)

  
!!$  write (*,*) 'IA', IA(:)
!!$  write (*,*) 'JA', JA(:)
!!$  write (*,*) 'AA', AA(:)
!!$  write (*,*) 'UPTR', UPTR(:)
!!$
!!$  pause
!!$
!!$  do i = 1, N
!!$    
!!$     k1 = ia(i)
!!$     k2 = ia(i+1)-1
!!$
!!$     do j=1, N
!!$        if (JA(k1) == j) then
!!$           write (*,*) AA(JA(k1))
!!$           k1 = k1 + 1
!!$        else
!!$           write (*,*) 0.0
!!$        end if
!!$     end do    
!!$
!!$     pause
!!$  end do
  
  
  call Timers_stop("diff_computeblkAMat")    
  
end subroutine diff_computeblkAMat
