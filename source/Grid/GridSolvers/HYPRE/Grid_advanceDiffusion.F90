!!****if* source/Grid/GridSolvers/HYPRE/Grid_advanceDiffusion
!!
!! NAME
!!
!!  Grid_advanceDiffusion
!!
!! SYNOPSIS
!!
!!  call Grid_advanceDiffusion (integer(IN)  :: iVar,
!!                              integer(IN)         :: iSrc,
!!                              integer(IN)  :: iFactorB,
!!                              integer(IN)  :: iFactorA,
!!                              integer(IN)  :: bcTypes(6),
!!                              real(IN)     :: bcValues(2,6),
!!                              real(IN)     :: dt,
!!                              real(IN)     :: chi,
!!                              real(IN)     :: scaleFact,
!!                              real(IN)     :: theta,
!!                              logical(IN)  :: solnIsDelta,
!!                     OPTIONAL,integer(IN)  :: iFactorC,
!!                     OPTIONAL,integer(IN)  :: iFactorD,
!!                     OPTIONAL,integer(IN)  :: pass)
!!
!!  DESCRIPTION 
!!
!!      This routine advances a generalized diffusion operator of the form
!!
!!         A*(df/dt) + C*f = div(B*grad(f)) + D ,
!!
!!      where
!!         f = f(x,t) is the  Variable to be diffused (x=1D..3D position);
!!
!!         A,B,C,D are optional given scalar factors/terms that may depend
!!         on position; they are either physcially constant in time, or at
!!         least considered time-independent for the purpose of the operation
!!         implemented here (typically by computing their values from the
!!         solution state reached by the previous time step).
!!
!!      Presently it is used to do heat conduction and multigroup diffusion.
!!
!! ARGUMENTS
!!   iVar           : Variable on which the diffusion operation is performed (e.g., TEMP_VAR)
!!   iFactorA       :| Are factors in the equation with spatial variation.
!!   iFactorB       :| Factor C,D are optional and are generally used
!!   iFactorC       :| to represent emission/absorption in MGD.
!!   iFactorD       :| iFactorA is not needed (and thus set to constant 1.0) for
!!                   | radiation diffusion.
!!                   | The factor B (whose role is that of diffusion coefficient) is passed to
!!                   | the solver implementation in cell-centered form in UNK; the solver will
!!                   | use averaging of adjacent cells (or perhaps some other mechanism) to
!!                   | derive the required face-cented diffusion coefficent values.
!!   bcTypes        : Presently OUTFLOW, VACUUM, DIRICHLET are supported, with additional
!!                    limited support for OUTSTREAM.
!!   bcValues       : Values of iVar,iFactorB on boundary (DIRICHLET).                        
!!   dt             : The time step.
!!   scaleFact      : Factor by which the end solution is scaled (not used).
!!   chi            : useful for constant diffusion problems (not used).
!!   theta          : varies scheme (0-> Explicit, 1-> backward euler, 0.5 -> Crank Nicholson
!!   pass           : Ignored in unsplit solver.
!!                    pass=1 order of directional sweep X-Y-Z, 
!!                    pass=2 order of directional sweep Z-Y-X.
!!   iSrc           : Ignored.
!!   solnIsDelta    : Is the solution only a delta that the caller has to apply to the
!!                    temperature, rather than temperature itself (ignored).
!!
!! SIDE EFFECTS
!!
!!
!! NOTES
!!  This implementation
!!    * supports: 
!!           1D, 2D, 3D Cartesian PARAMESH (with local refinement)
!!           1D, 2D, 3D Cartesian in UG (3D untested).
!!           1D Spherical in PARAMESH/UG
!!           1D, 2D Cylindrical in PARAMESH/UG (R-Z)
!!
!!    *  uses HYPRE library to solve AX = B
!!
!!
!!***

!!REORDER(4): solnVec

subroutine Grid_advanceDiffusion (iVar, iSrc, iFactorB, iFactorA, bcTypes, bcValues, &
     dt, chi, scaleFact,theta, solnIsDelta,iFactorC, iFactorD, pass)
  
  use Grid_data,        ONLY : gr_meshMe, gr_meshcomm
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_interface,     ONLY : gr_hypreCreateMatrix, gr_hypreComputeB,    &
                               gr_hypreGridStatus
  use gr_hypreLocalInterface, ONLY: gr_hypreExchangeFacB
  use Grid_interface,   ONLY : Grid_fillGuardCells, Grid_getListOfBlocks, &
                               Grid_getBlkPtr, Grid_releaseBlkPtr,        &
                               Grid_getBlkIndexLimits, Grid_getBlkData,   &
                               Grid_getBlkRefineLevel
  
  use gr_hypreData,   ONLY   : gr_hypreLower, gr_hypreUpper, &
                               gr_hypreMatA, gr_hypreVecB, gr_hypreRefineMIN, &
                               gr_hyprePrintSolveInfo, &
                               gr_hypreFloor, gr_hypreRelTol, &
                               gr_hypreAbsTol, gr_hypreSolverAbsTolEff, &
                               gr_hypreSolverAutoAbsTolFact

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Flash_mpi.h"
  
  integer, intent(IN) :: iVar
  integer, intent(IN) :: iSrc
  integer, intent(IN) :: iFactorB
  integer, intent(IN) :: iFactorA
  real, intent(IN)    :: dt 
  real, intent(IN)    :: chi
  real, intent(IN)    :: scaleFact
  real, intent(IN)    :: theta
  logical, intent(IN) :: solnIsDelta
  integer, dimension(6),  intent(IN) :: bcTypes
  real   , dimension(2,6),intent(IN) :: bcValues
  integer, intent(IN), OPTIONAL :: pass
  integer, intent(IN), OPTIONAL :: iFactorC
  integer, intent(IN), OPTIONAL :: iFactorD   
  
  integer :: blockCount
  integer :: blockList(MAXBLOCKS)

  
!!! BEGIN TEMP

  integer :: mylevel, mypart, var, ii,i,j,k, ierr
  real, allocatable :: BoxVal(:),prevBoxVal(:)
  real, allocatable :: RHSVal(:)
  integer :: datasize(MDIM)
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  integer :: blockID, lb
  real, allocatable :: cellVolumes(:,:,:)
  integer, dimension(2,MDIM):: blkLimitsGC, blkLimits 
  
  character(len=32) :: matfile

!!! END TEMP

  real :: normAllDiagEst2, normsAllDiagEst(2:2)
  logical :: estimPrec
  real :: relevantNorm
  
  call Timers_start("Grid_advanceDiffusion")   
  
!!$  call Timers_start("Diffusion Barrier")
!!$  call MPI_Barrier (gr_meshcomm, ierr)
!!$  call Timers_stop ("Diffusion Barrier")
  
  !print*,theta; stop

  call Grid_getListOfBlocks(LEAF,blockList,blockCount)
  
  !!-----------------------------------------------------------------------
  !!     1.  Do we need to reset HYPRE grid ?, has the underlying AMR
  !!         mesh been modified ?, is this the first call to HYPRE ?
  !!         Only in AMR is this routine meaningful.
  !!-----------------------------------------------------------------------
  call gr_hypreGridStatus (blockCount, blockList)  
  
  !!-----------------------------------------------------------------------
  !!     2. Exchange iFactorB across processors. Needs to be done only in
  !!        PARAMESH / AMR, if UG the function will return without any action.
  !!-----------------------------------------------------------------------
  call gr_hypreExchangeFacB (iFactorB, blockCount, blockList)  
  
  !!-----------------------------------------------------------------------
  !!     3. Set initial guess for solver typically 
  !!        iVar at previous time step is used.
  !!-----------------------------------------------------------------------
  call gr_hypreSetIniGuess (iVar, blockCount, blockList)  
    
  
  !!-----------------------------------------------------------------------
  !!     4. Compute div B grad terms of a matrix.
  !!-----------------------------------------------------------------------
  if (theta /= 0.0) then !! Implicit
     !!     4a. Compute (part of) the actual A matrix in AX=B
     call gr_hypreCreateMatrix(iVar, iFactorB, iFactorA, bcTypes, bcValues, dt, theta,  &
       blockCount, blockList, .TRUE.)     
  else
     !!     4b. Assemble (part of) a matrix M such that B = MX, where B is the RHS of AX=B
     !!         and X is the initial guess, i.e., iVar from the previous time step.
     call gr_hypreCreateMatrix(iVar, iFactorB, iFactorA, bcTypes, bcValues, dt, theta, &
       blockCount, blockList, .FALSE.)
  end if
  
  !!-----------------------------------------------------------------------
  !!     5. Using matrix M compute B = MX; additionally if needed, iFactorD
  !!        is to be added to the matrix diagonal. Doing this step helps 
  !!        avoid doing a point to point communication to share X. This is 
  !!        now done by HYPRE implicitly.
  !!-----------------------------------------------------------------------
  call gr_hypreComputeB (blockCount, blockList, iVar, iFactorA, iFactorB, dt, theta, &
       bcTypes, bcValues, iFactorD)     
  
  !!-----------------------------------------------------------------------
  !!     4c. Compute part of the actual A matrix in AX=B if not done in 4a.
  !!-----------------------------------------------------------------------
  if (theta == 0.0) then
     !!     Actually, this should be equivalent to just zeroing out the matrix!
     call gr_hypreCreateMatrix(iVar, iFactorB, iFactorA, bcTypes, bcValues, dt, 0.0,  &
       blockCount, blockList, .TRUE.)
  end if
  
  
  !!-----------------------------------------------------------------------
  !!     6.  Finalize computation of the actual A matrix in AX=B, all terms.
  !!         - add term from Factor A to diagonal
  !!         - add term from Factor C to diagonal
  !!     6a. Finalize computation of the RHS vector B in AX=B, all terms.
  !!         - add term from Factor A
  !!         - (do not add term from Factor C - not currently included in RHS)
  !!         - add term from "Factor D"
  !!-----------------------------------------------------------------------

  mypart = 0  !! part iterator 
  var    = 0  !! var iterator.
  
  estimPrec = (gr_hypreSolverAutoAbsTolFact .NE. 0.0)

  normAllDiagEst2 = 0.0
111 format(A,1x,ES20.13)

  do lb = 1, blockCount 
     
     blockID = blockList(lb)
     
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)    
     call Grid_getBlkPtr(blockID, solnVec)
     call Grid_getBlkRefineLevel(blockID,mylevel)

     mypart = mylevel - gr_hypreRefineMIN
     
     datasize(1:MDIM) = blkLimits(HIGH,1:MDIM)-blkLimits(LOW,1:MDIM)+1         
     
     allocate(cellVolumes(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)))
     
     call Grid_getBlkData(blockID, CELL_VOLUME, 0, EXTERIOR,          &
          blkLimits(LOW,:), cellVolumes,datasize)      
     
     allocate(BoxVal(product(dataSize(1:NDIM))))
     allocate(prevBoxVal(product(dataSize(1:NDIM))))
     allocate(RHSVal(product(dataSize(1:NDIM))))     
     
     if (estimPrec) then
        call HYPRE_SStructMatrixGetBoxValues(gr_hypreMatA, mypart, gr_hypreLower(lb, 1:NDIM), & 
             gr_hypreUpper(lb,1:NDIM), var, 1, (/0/), prevBoxVal(:), ierr)
        if (ierr .NE. 0) then
           call Driver_abortFlash('Got nonzero ierr from HYPRE_SStructMatrixGetBoxValues!')
        end if
     end if

     do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)      
        do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
           do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
              
              ii = (k - blkLimits(LOW,KAXIS)  + 1)                             +  &
                   (j - blkLimits(LOW,JAXIS))*dataSize(KAXIS)                  +  &
                   (i - blkLimits(LOW,IAXIS))*dataSize(KAXIS)*dataSize(JAXIS)  
              
              BoxVal(ii) = solnVec(iFactorA,i,j,k)*cellVolumes(i,j,k)
              
              RHSVal(ii) = BoxVal(ii)*solnVec(iVar,i,j,k)
              
              if (present(iFactorC)) then
                 BoxVal(ii) = BoxVal(ii) + dt*solnVec(iFactorC,i,j,k)*cellVolumes(i,j,k)         
                 if (estimPrec) then
                    if ((BoxVal(ii) .LE. 0.0)) then
                       print*,'BAD DIAG VALUE!!!', blockid,ii, &
                            BoxVal(ii),BoxVal(ii)- dt*solnVec(iFactorC,i,j,k)*cellVolumes(i,j,k)
                    end if
                 end if
              end if
              
              if (present(iFactorD)) then 
                 RHSVal(ii) =  RHSVal(ii) + SolnVec(iFactorD,i,j,k)*cellVolumes(i,j,k)*dt
              end if
              if (estimPrec) &
                   normAllDiagEst2 = normAllDiagEst2 + gr_hypreFloor**2 / (BoxVal(ii)+prevBoxVal(ii))
                            
           end do
        end do
     end do
     
     
     call HYPRE_SStructVectorAddToBoxValu(gr_hypreVecB, mypart, gr_hypreLower(lb, 1:NDIM), &
          gr_hypreUpper(lb,1:NDIM), var, RHSVal(:), ierr)
     
     
     call HYPRE_SStructMatrixAddToBoxValu(gr_hypreMatA, mypart, gr_hypreLower(lb, 1:NDIM), & 
          gr_hypreUpper(lb,1:NDIM), var, 1, (/0/), BoxVal(:), ierr)     
     
     
     call Grid_releaseBlkPtr(blockID, solnVec)
     
     deallocate (BoxVal)
     deallocate (prevBoxVal)
     deallocate (RHSVal)
     deallocate(cellVolumes)
     
  end do  
  
  
!!$  matfile = 'ex12f.out'
!!$  matfile(10:10) = char(0)
!!$  call HYPRE_SStructMatrixPrint(matfile, gr_hypreMatA, 0, ierr)
!!$  pause

  if (estimPrec) then
     call MPI_ALLREDUCE(normAllDiagEst2,normsAllDiagEst,1,FLASH_REAL,MPI_SUM,&
          gr_meshComm,ierr)
  end if

  if (estimPrec) then
     normsAllDiagEst(2) = sqrt(normsAllDiagEst(2))
     if (gr_meshMe==MASTER_PE .AND. .FALSE.) then
        print 111,"AllDiag-C norm estimates:   ", normsAllDiagEst(2)
     end if
  end if
  

  if (estimPrec) then
     normsAllDiagEst(2) =      normsAllDiagEst(2) * gr_hypreSolverAutoAbsTolFact
     relevantNorm = normsAllDiagEst(2)
     if (gr_meshMe==MASTER_PE .AND. gr_hyprePrintSolveInfo) then
        print 111,"proposed abs tol (norm ALL):", normsAllDiagEst(2)
     end if

     gr_hypreSolverAbsTolEff = relevantNorm

     if (gr_hypreAbsTol > 0.0) gr_hypreSolverAbsTolEff = max(gr_hypreSolverAbsTolEff, gr_hypreAbsTol)
  else
     if (gr_hypreAbsTol > 0.0) gr_hypreSolverAbsTolEff = gr_hypreAbsTol
  end if
  
  !!-----------------------------------------------------------------------
  !!     7. Solve AX = B
  !!-----------------------------------------------------------------------
  call gr_hypreSolve ()
  
  !!-----------------------------------------------------------------------
  !!     8. Update unk variable using X.
  !!-----------------------------------------------------------------------
  call gr_hypreUpdateSoln (iVar, blockCount, blockList)  
  
  
  call Timers_stop("Grid_advanceDiffusion") 
  
  return
  
end subroutine Grid_advanceDiffusion
