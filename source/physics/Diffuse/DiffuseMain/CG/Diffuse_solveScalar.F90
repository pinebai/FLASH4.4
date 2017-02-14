!!****if* source/physics/Diffuse/DiffuseMain/CG/Diffuse_solveScalar
!!
!!  NAME 
!!
!!  Diffuse_solveScalar
!!
!!  SYNOPSIS
!!
!!  call Diffuse_solveScalar (integer, intent(IN) :: iVar,
!!                            integer, intent(IN) :: iFactorB,
!!                            integer, intent(IN) :: iFactorA,
!!                            integer, intent(IN) :: bcTypes(6),
!!                            real,    intent(IN) :: bcValues(2,6),
!!                            real,    intent(IN) :: dt,
!!                            real,    intent(IN) :: scaleFact,
!!                            real,    intent(IN) :: chi,
!!                            real,    intent(IN) :: theta,
!!                            integer, OPTIONAL, intent(IN) :: pass,
!!                            integer, intent(IN) :: blockCount,
!!                            integer,dimension(blockCount),intent(IN) :: blockList,
!!                            integer, intent(IN), OPTIONAL :: iFactorC,
!!                            integer, intent(IN), OPTIONAL :: iFactorD)
!!
!!
!!  DESCRIPTION 
!!
!!      This routine advances the diffusion operator of the form,
!!      A*(df/dt) + C*f = div(B*grad(f)) + D
!!      f -> Variable to be diffused.
!!      C,D are optional factors.
!!  
!!      Presently it is used to do conduction and multigroup diffusion.
!!
!! ARGUMENTS
!! 
!!   iVar           : Variable on which the diffusion operatorion is performed (e.g TEMP_VAR)
!!   iFactorA       :| Are factors in the equation with spatial variation.
!!   iFactorB       :| Factor C,D are optional and are generally used
!!   iFactorC       :| to represent emission/absorption in MGD.
!!   iFactorD       :| iFactorA is needed only for conduction.
!!   bcTypes        : Presently OUTFLOW, VACUUM is supported, DIRICHLET is untested.
!!   bcValues       : Values of iVar,iFactorB on boundary (DIRICHLET).                        
!!   dt             : The time step.
!!   scaleFact      : Factor by which the end solution is scaled (not used).
!!   chi            : useful for constant diffusion problems (not used).
!!   theta          : varies scheme (0-> Explicit, 1-> backward euler, 0.5 -> Crank Nicholson
!!   pass           : Ignored in unsplit solver.
!!                    pass=1 order of directional sweep X-Y-Z, 
!!                    pass=2 order of directional sweep Z-Y-X.
!!   blockCount     : The number of blocks in the list.   
!!   blockList      : The list of blocks on which the solution must be updated.                    
!!
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES
!!
!!  Stub implementation.              
!!
!!  
!!
!!***

!!REORDER(4): solnVec

#include "Flash.h"

subroutine Diffuse_solveScalar (iVar, iFactorB, iFactorA, bcTypes, bcValues, &
     dt, scaleFact, chi, theta, pass, blockCount, blockList, iFactorC, iFactorD)
  
  use diff_interface, ONLY : diff_computeAX, diff_computeblkAMat
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
                             Grid_getBlkIndexLimits, Grid_fillGuardCells, &
                             Grid_getBlkData
  use Diffuse_data,   ONLY : diff_meshMe, diff_geometry, diff_meshComm
  
  
  implicit none
  
#include "constants.h"  
#include "Flash_mpi.h"
  
  integer, intent(IN):: iVar
  integer, intent(IN):: iFactorB
  integer, intent(IN):: iFactorA
  integer, intent(IN):: bcTypes(6)
  integer, intent(IN):: blockCount
  integer, dimension(blockCount),intent(IN):: blockList
  real, intent(IN):: bcValues(2,6)
  real, intent(IN):: dt
  real, intent(IN):: scaleFact
  real, intent(IN):: chi
  real, intent(IN):: theta
  integer, OPTIONAL,intent(IN):: pass
  integer, OPTIONAL,intent(IN):: iFactorC
  integer, OPTIONAL,intent(IN):: iFactorD

  !====================================================================================================
  real, allocatable, dimension(:,:,:,:):: AP, R, Z, B
  logical :: mask(NUNK_VARS)
  integer, dimension(MDIM) :: maxSize    
  integer, dimension(2,MDIM):: blkLimitsGC, blkLimits
  integer :: lb, blockID
  real, POINTER, DIMENSION(:,:,:,:) :: solnData
  integer :: i,j,k
  logical :: converged
  real, dimension(blockCount) :: RdotZ, APdotP,newRdotZ
  real :: GRdotZ, GAPdotP, GnewRdotZ
  real :: initialRes, alpha, beta
  integer :: ierr, iter

  
  !! for CSR A Matrix (ILU(0) AS PC)
  integer :: N,NZ
  integer, allocatable,dimension(:)   :: IA, JA, UPTR
  real,    allocatable,dimension(:,:) :: AA

  logical :: diff_usePrecon = .true.  
    
  call Grid_fillGuardCells(CENTER,ALLDIR) 
  
  mask = .false.
  mask(KVEC_VAR) = .true.  
  
  maxSize = 0
  do lb = 1, blockCount           
     call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)                      
     maxSize(:) = MAX (blkLimitsGC(HIGH,:),maxSize(:))    
  end do
  
  !! RESIDUE.
  allocate(R(blockCount,maxSize(IAXIS),maxSize(JAXIS),maxSize(KAXIS)))
  
  !! MZ = R.
  allocate(Z(blockCount,maxSize(IAXIS),maxSize(JAXIS),maxSize(KAXIS))) 
  
  !! MATRIX_VECTOR_PRODUCT(A,P).
  allocate(AP(blockCount,maxSize(IAXIS),maxSize(JAXIS),maxSize(KAXIS))) 
  
  !! RHS VECTOR.
  allocate(B(blockCount,maxSize(IAXIS),maxSize(JAXIS),maxSize(KAXIS))) 
  
  !! INITIAL GUESS FOR SOLVER (DOES FLUX CONSERVATION)
  call diff_computeFluxes (blockCount,blockList,iVar,iFactorB)  
  
  !! COMPUTE RHS
  do lb = 1, blockCount         
     blockID = blockList(lb)
     
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)                      
     
     call diff_computeAX (blockID, blkLimits, blkLimitsGC, iVar, & 
          iFactorA, dt, (theta-1.0), B(lb,:,:,:),iFactorC,iFactorD)
     
  end do


  !! PRECONDITIONER
  maxSize(IAXIS) = maxSize(IAXIS) - 2*NGUARD
  maxSize(JAXIS) = maxSize(JAXIS) - 2*NGUARD*K2D
  maxSize(KAXIS) = maxSize(KAXIS) - 2*NGUARD*K3D
  
  N = PRODUCT(maxSize)
  
  if (diff_usePrecon) then    
     
     if (NDIM == 1) then
        NZ = 3*(maxSize(IAXIS)-2) + 2*2
     else if (NDIM == 2) then
        NZ = 5*(maxSize(IAXIS)-2)*(maxSize(JAXIS)-2) + & 
             4*(2*(maxSize(JAXIS)-2)+2*(maxSize(IAXIS)-2)) + 3*4  
     else        
        NZ = 7*(maxSize(IAXIS)-2)*(maxSize(JAXIS)-2)*(maxSize(KAXIS)-2) + &
             5*4*(maxSize(IAXIS)+maxSize(JAXIS)+maxSize(KAXIS)-6) + &
             6*2*((maxSize(IAXIS)-2)*(maxSize(JAXIS)-2) + (maxSize(IAXIS)-2)*(maxSize(KAXIS)-2) & 
             + (maxSize(JAXIS)-2)*(maxSize(KAXIS)-2))+&
             4*8        
     end if
     
     allocate (JA(NZ))
     allocate (IA(N+1))
     allocate (UPTR(N))    
     allocate (AA(NZ,blockCount))    
     
     do lb = 1, blockCount         
        
        blockID = blockList(lb)
        
        call diff_computeblkAMat (blockID, NZ, N, dt, theta, AA(:,lb), & 
             JA, IA, UPTR, iFactorA, iFactorB, iFactorC) 
        
        call diff_computeILU(AA(:,lb), JA, IA, UPTR, NZ, N)  !! ILU is stored in AA.        
        
     end do
     
  end if
  
  AP = 0.0

  do lb = 1, blockCount         
     
     blockID = blockList(lb)

     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)    
     
     call Grid_getBlkPtr(blockID, solnData)           
     
     !!COMPUTE AX
     call diff_computeAX (blockID, blkLimits, blkLimitsGC, iVar, iFactorA, & 
          dt, theta, AP(lb,:,:,:),iFactorC,iFactorD)    
     
     !! COMPUTE RESIDUE
     R(lb,:,:,:) = B(lb,:,:,:)- AP(lb,:,:,:)
     
     if (diff_usePrecon) then 
        call diff_applyPC (Z(lb,:,:,:), R(lb,:,:,:), NZ, N, AA(:,lb), & 
             JA, IA, UPTR, blkLimits, blkLimitsGC)     
     else
        !! Z /= R if a PC is used.
        Z(lb,:,:,:) = R(lb,:,:,:) 
     endif
     
     ! P0 = Z0, FIRST KRYLOV VECTOR.
     solnData(KVEC_VAR,:,:,:) = Z(lb,:,:,:)    
     
     call diff_dotproduct(R(lb,:,:,:), Z(lb,:,:,:), RdotZ(lb), & 
          blkLimits, blkLimitsGC)     ! rj dot zj  
     
     call Grid_releaseBlkPtr(blockID, solnData)
     
  end do
  
  converged = .false.  
  
  call mpi_allreduce (sum(RdotZ),GRdotZ, 1, FLASH_REAL, MPI_SUM, diff_meshComm, ierr)   
  
  initialRes = GRdotZ  
  
  iter = 0
  
  if (sqrt(initialRes) > 1.0e-12) then     
     
     do while (.not. converged)              
        
        call Grid_fillGuardCells(CENTER,ALLDIR,masksize=NUNK_VARS, & 
             mask=mask, selectBlockType=LEAF)
        
        APdotP   = 0.
        GAPdotP  = 0.
        AP       = 0.0
        
        call diff_computeFluxes(blockCount,blockList,KVEC_VAR,iFactorB)  
        
        do lb = 1, blockCount                   

           blockID = blockList(lb)                     

           call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)               

           call diff_computeAX (blockID,blkLimits,blkLimitsGC,KVEC_VAR, &
                iFactorA,dt,theta, AP(lb,:,:,:),iFactorC,iFactorD)
           
           call Grid_getBlkPtr(blockID, solnData)                       
           
           call diff_dotproduct(AP(lb,:,:,:), solnData(KVEC_VAR,:,:,:), & 
                APdotP(lb),blkLimits, blkLimitsGC)            
           
           call Grid_releaseBlkPtr(blockID, solnData)                   
           
        end do
        
        call mpi_allreduce (sum(APdotP), GAPdotP, 1, FLASH_REAL, & 
             MPI_SUM, diff_meshComm, ierr)
        
        alpha = GRdotZ/GAPdotP            
        
        do lb = 1, blockCount      
           
           blockID = blockList(lb)
           
           call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)           
           
           call Grid_getBlkPtr(blockID, solnData)           
           
           solnData(iVar,:,:,:) =  solnData(iVar,:,:,:)  + alpha*solnData(KVEC_VAR,:,:,:)                  
           
           R(lb,:,:,:) = R(lb,:,:,:) - alpha*AP(lb,:,:,:)
           
           if (diff_usePrecon) then                  
              call diff_applyPC (Z(lb,:,:,:), R(lb,:,:,:), NZ, N, &
                   AA(:,lb), JA, IA, UPTR, blkLimits, blkLimitsGC)
           else
              
              Z(lb,:,:,:) = R(lb,:,:,:)
              
           end if
           
           !!Inner product of newly computed residue.
           call diff_dotproduct(R(lb,:,:,:),Z(lb,:,:,:),newRdotZ(lb),blkLimits, blkLimitsGC)       
           
           call Grid_releaseBlkPtr(blockID, solnData)
           
        end do
        
        call mpi_allreduce (sum(newRdotZ), GnewRdotZ, 1, FLASH_REAL, &
             MPI_SUM, diff_meshComm, ierr)   
        
        beta = GnewRdotZ / GRdotZ   
        
       do lb = 1, blockCount       
          
          blockID = blocklist(lb)
          
          call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)           
          
          call Grid_getBlkPtr(blockID,solnData)                                     
          
          solnData(KVEC_VAR,:,:,:) = Z(lb,:,:,:) + beta*solnData(KVEC_VAR,:,:,:)          
          
          call Grid_releaseBlkPtr(blockID,solnData)    
          
       end do
       
       GRdotZ = GnewRdotZ        
       
       iter = iter + 1              
       
       !! check for convergence, diverenge or Max iterations.
       if (sqrt(GRdotZ) <= MAX(1.0e-12*sqrt(initialRes), 1.0e-50) .or. iter > 1000) then
          converged = .true.  
       else
          if (sqrt(GRdotZ) >= 1.0e3*sqrt(initialRes)) then                          
             call Driver_abortFlash("Conjugate Gradient residue is diverging!") 
          endif
       end if
       
    end do
    
 end if

 
 !!Check for Max iterations, abort if more then set number of iterations.
 if ( diff_meshMe == MASTER_PE ) then
    if (iter > 1000) then
       call Driver_abortFlash("Conjugate Gradient: No convergence after 1000 iterations, try higher value!") 
    end if
    
    !write (*,*) iter, sqrt(GRdotZ/initialRes), blockCount
    
 end if
 
 !! Flooring, avoids small (<< 0.0) negative floating point numbers.
 do lb = 1, blockCount
    call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)
    call Grid_getBlkPtr(blocklist(lb), solnData)    
    solnData(iVar,:,:,:) = max(solnData(iVar,:,:,:), 1.0E-10) 
    call Grid_releaseBlkPtr(blocklist(lb), solnData)
 end do
 
 if (diff_usePrecon) then
    deallocate (JA)
    deallocate (IA)
    deallocate (UPTR)    
    deallocate (AA)  
 end if
 
 deallocate(B)
 deallocate(R)
 deallocate(Z)
 deallocate(AP)  
 
 return
  
end subroutine Diffuse_solveScalar




subroutine diff_dotproduct (vec1, vec2, dotprod, blkLimits, blkLimitsGC)   
  
  implicit none

#include "constants.h"   
  
  integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits, blkLimitsGC
  
  real, intent(IN)    :: vec1(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))  
  real, intent(IN)    :: vec2(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))
  
  real, intent (OUT)  :: dotprod
  
  integer :: i, j, k
  
  dotprod = 0.    
  
  do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
     do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
           dotprod = dotprod + vec1(i,j,k)*vec2(i,j,k)
        end do
     end do
  end do


end subroutine diff_dotproduct
