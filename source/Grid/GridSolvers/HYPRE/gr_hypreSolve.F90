!!****if* source/Grid/GridSolvers/HYPRE/gr_hypreSolve
!!
!!  NAME 
!!
!!  gr_hypreSolve
!!
!!  SYNOPSIS
!!
!!  call gr_hypreSolve ()
!!
!!  DESCRIPTION 
!!      This routine solves AX=B using HYPRE
!!
!! ARGUMENTS
!!
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES
!!
!!      Requires HYPRE library.
!!
!!***


subroutine gr_hypreSolve()
  
  use gr_hypreData,     ONLY : gr_hypreMatA, gr_hypreVecB, gr_hypreVecX, &
                               gr_hypreSolver, gr_hypreSolverType, &
                               gr_hyprePrintSolveInfo,gr_hypreMaxIter, &
                               gr_hypreRelTol, gr_hypreSolverAbsTolEff
  use gr_solversData,   ONLY : dbgContextSolvers => gr_solversDbgContext
  use Timers_interface, ONLY : Timers_start, Timers_stop  
  use Driver_interface, ONLY : Driver_getNStep
  use RadTrans_interface,ONLY: RadTrans_dbgContext_t,RadTrans_getDbgContextPtr
  use Grid_data,        ONLY : gr_meshMe
  use gr_hypreLocalInterface, ONLY : hypre_pcggetconverged, hypre_describeerror
  
  implicit none
  
#include "Flash.h"  
#include "constants.h"
#include "HYPREf.h"   

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif
#ifdef DEBUG_GRID
#define DEBUG_SOLVERS
#endif
#ifdef DEBUG_SOLVERS
#define DEBUG_MATRIX
#define DEBUG_VECTOR
#endif

  integer, parameter :: range=SELECTED_INT_KIND(16)

  integer (KIND=range)::  parA
  integer (KIND=range)::  parb
  integer (KIND=range)::  parx 
  integer :: ierr
  integer :: nstep
  integer :: hypreError
  integer :: iConverged
  logical :: converged, convergedAndOk

  integer :: num_iterations
  real    :: final_res_norm, initial_rhs_norm
  integer :: component, group
  type(RadTrans_dbgContext_t),pointer :: dbgContextRt

  integer, parameter :: debugMatMinStep=218, debugMatMaxStep=219
  integer, parameter :: debugMatMinGroup=5, debugMatMaxGroup=6
  character(len=22) :: matfile
  character(len=7)  :: nstepStr
  character(len=5)  :: groupStr
  

  call Timers_start("gr_hypreSolve")     
  
  call HYPRE_SStructVectorAssemble(gr_hypreVecB, ierr)  

  component        = -1
  group            = -1
  num_iterations   = 0  
  final_res_norm   = 0.0
  converged        = .TRUE.     !for other solver than HYPRE_PCG
  convergedAndOk   = .FALSE.
  iConverged       = -1

  !! Rob Falgout recommended on 2013-01-09 to do the solver setup
  !! each time before the solver is called.
  call gr_hypreSetupSolver()

  if (gr_hypreSolverType == HYPRE_SPLIT) then     
     call HYPRE_SStructSplitSetup(gr_hypreSolver, gr_hypreMatA, gr_hypreVecB, gr_hypreVecX, ierr)
     call HYPRE_SStructSplitSolve(gr_hypreSolver, gr_hypreMatA, gr_hypreVecB, gr_hypreVecX, ierr)
  else  

     component = dbgContextSolvers%component
     group = dbgContextSolvers%group
#if defined(DEBUG_MATRIX) || defined(DEBUG_VECTOR)
     call Driver_getNStep(nstep)
     write(nstepStr,'(I6.2)') nstep
     write(groupStr,'(I4.3)') dbgContextSolvers%group
     nstepStr = 'N' // adjustl(trim(nstepStr))
     groupStr = 'G' // adjustl(trim(groupStr))
111  format(A4,A,A5,A)
     if (nstep .GE. debugMatMinStep .AND. nstep .LE. debugMatMaxStep) then
        if (group .GE. debugMatMinGroup .AND. group .LE. debugMatMaxGroup) then

#ifdef DEBUG_MATRIX
           write(matfile,FMT=111) 'amat',trim(groupStr),'b4slv',nstepStr
           matfile(index(matfile,' '):) = char(0)
           call HYPRE_SStructMatrixPrint(trim(matfile), gr_hypreMatA, 1, ierr)
#endif

#ifdef DEBUG_VECTOR
           write(matfile,FMT=111) 'bvec',trim(groupStr),'b4slv',nstepStr
           matfile(index(matfile,' '):) = char(0)
           call HYPRE_SStructVectorPrint(trim(matfile), gr_hypreVecB, 0, ierr)
#endif

        end if
     end if
#endif

     !!-----------------------------------------------------------------------
     !! Because we are using a PARCSR solver, we need to get the object
     !! of the matrix and vectors to pass in to the ParCSR solvers.
     !!-----------------------------------------------------------------------  
     call HYPRE_SStructMatrixGetObject(gr_hypreMatA, parA, ierr)
     call HYPRE_SStructVectorGetObject(gr_hypreVecB, parb, ierr)
     call HYPRE_SStructVectorGetObject(gr_hypreVecX, parx, ierr)
     
        call Hypre_sstructinnerprod(gr_hypreVecB,gr_hypreVecB,initial_rhs_norm,ierr)

     if (gr_hypreSolverType == HYPRE_PCG) then            
        
!!$        call Hypre_sstructinnerprod(gr_hypreVecB,gr_hypreVecB,initial_rhs_norm,ierr)
#ifdef DEBUG_SOLVERS
        if (gr_meshMe == 0) then
           print*,'RHS_initial_norm**2 is',initial_rhs_norm
        end if
#endif
        if (ierr .NE. 0) then
           print*,'Error from Hypre_sstructinnerprod, continuing anyway on',gr_meshMe
        end if

        if (initial_rhs_norm .NE. 0.0) then
           call HYPRE_ParCSRPCGSetATol(gr_hypreSolver, gr_hypreSolverAbsTolEff, ierr)
           call Timers_start("HYPRE_ParCSRPCGSetup")    
           call HYPRE_ParCSRPCGSetup(gr_hypreSolver, parA, parb, parx, ierr) 
!!$           call HYPRE_ParCSRPCGSetPrintLevel(gr_hypreSolver, 2, ierr) 
           call Timers_stop("HYPRE_ParCSRPCGSetup")   

           call Timers_start("HYPRE_ParCSRPCGSolve")    
           call HYPRE_ParCSRPCGSolve(gr_hypreSolver, parA, parb, parx, ierr)     
           call Timers_stop("HYPRE_ParCSRPCGSolve")        
!!$           call HYPRE_SStructVectorGather(gr_hypreVecX, ierr)  
!!$           call hypre_sstructvectorsetconstantv(gr_hypreVecX,100.0,ierr)
        else
!!$           call hypre_sstructvectorsetconstantv(gr_hypreVecX,10.0,ierr)
        end if
        
     else  if (gr_hypreSolverType == HYPRE_BICGSTAB) then
        
        call Timers_start("HYPRE_ParCSRBiCGSTABSetup")    
        call HYPRE_ParCSRBiCGSTABSetup(gr_hypreSolver, parA, parb, parx, ierr)      
        call Timers_stop ("HYPRE_ParCSRBiCGSTABSetup")    
        
        call Timers_start("HYPRE_ParCSRBiCGSTABSolve")    
        call HYPRE_ParCSRBiCGSTABSolve(gr_hypreSolver, parA, parb, parx, ierr)     
        call Timers_stop("HYPRE_ParCSRBiCGSTABSolve")    
        
     else if(gr_hypreSolverType == HYPRE_AMG) then         
        
        call Timers_start("HYPRE_BoomerAMGSetup")
        call HYPRE_BoomerAMGSetup(gr_hypreSolver, parA, parb, parx, ierr)
        call Timers_stop("HYPRE_BoomerAMGSetup")
        
        call Timers_start("HYPRE_BoomerAMGSolve")
        call HYPRE_BoomerAMGSolve(gr_hypreSolver, parA, parb, parx, ierr)
        call Timers_stop("HYPRE_BoomerAMGSolve")                            
        
     else if(gr_hypreSolverType == HYPRE_GMRES) then   
        
        call Timers_start("HYPRE_ParCSRGMRESSetup")    
        call HYPRE_ParCSRGMRESSetup(gr_hypreSolver, parA, parb, parx, ierr)      
        call Timers_stop("HYPRE_ParCSRGMRESSetup")    
        
        call Timers_start("HYPRE_ParCSRGMRESSolve")    
        call HYPRE_ParCSRGMRESSolve(gr_hypreSolver, parA, parb, parx, ierr)
        call Timers_stop("HYPRE_ParCSRGMRESSolve")    
        
     else if(gr_hypreSolverType == HYPRE_HYBRID) then           
        
        call Timers_start("HYPRE_ParCSRHybridsetup")    
        call HYPRE_ParCSRHybridsetup(gr_hypreSolver,parA, parb, parx, ierr)
        call Timers_stop("HYPRE_ParCSRHybridsetup")    
        
        call Timers_start("HYPRE_ParCSRHybridSolve")    
        call HYPRE_ParCSRHybridSolve(gr_hypreSolver, parA, parb, parx, ierr)
        call Timers_stop("HYPRE_ParCSRHybridSolve")    
     end if
     
  end if

  if (gr_meshMe == 0 .OR. gr_hypreSolverType .NE. HYPRE_SPLIT) then
     
     select case (gr_hypreSolverType)
        
     case (HYPRE_PCG)
        if (initial_rhs_norm .NE. 0.0) then
           call HYPRE_ParCSRPCGGetNumIterations(gr_hypreSolver, num_iterations, ierr)
           call HYPRE_ParCSRPCGGetFinalRelative(gr_hypreSolver, final_res_norm, ierr)
!!$        print*, "HYPRE PCG Num Iterations = ", num_iterations 
!!$        print*, "HYPRE PCG Relative Residual Norm = ", final_res_norm
           call hypre_pcggetconverged(gr_hypreSolver, iConverged, ierr)
           converged = (iConverged .NE. 0)
        else
           num_iterations = 0
           final_res_norm = 0.0
           converged = .TRUE.
           iConverged = -2
        end if
        call parseSolveStatus
        
     case (HYPRE_AMG)        
        call HYPRE_BoomerAMGGetNumIterations(gr_hypreSolver, num_iterations, ierr)
        call HYPRE_BoomerAMGGetFinalReltvRes(gr_hypreSolver, final_res_norm, ierr)
        call parseSolveStatus
!!$        print*, "HYPRE AMG Num Iterations = ", num_iterations 
!!$        print*, "HYPRE AMG Relative Residual Norm = ", final_res_norm        
        
     case (HYPRE_BICGSTAB)                
        call HYPRE_ParCSRBICGSTABGetNumIter(gr_hypreSolver, num_iterations, ierr)
        call HYPRE_ParCSRBICGSTABGetFinalRel(gr_hypreSolver, final_res_norm, ierr)        
        call parseSolveStatus
!!$        print*, "HYPRE BiCGSTAB Num Iterations = ", num_iterations 
!!$        print*, "HYPRE BiCGSTAB Relative Residual Norm = ", final_res_norm
        
     case (HYPRE_GMRES)                
        call HYPRE_ParCSRGMRESGetNumIteratio(gr_hypreSolver, num_iterations, ierr)
        call HYPRE_ParCSRGMRESGetFinalRelati(gr_hypreSolver, final_res_norm, ierr)        
        call parseSolveStatus
!!$        print*, "HYPRE GMRES Num Iterations = ", num_iterations 
!!$        print*, "HYPRE GMRES Relative Residual Norm = ", final_res_norm
        
     case (HYPRE_HYBRID)           
        
        call HYPRE_ParCSRHybridGetNumIterati(gr_hypreSolver, num_iterations, ierr)
!!$        print*, "HYPRE HYBRID Num Iterations      = ", num_iterations 
        
        call HYPRE_ParCSRHybridGetPCGNumIter (gr_hypreSolver, num_iterations, ierr)
!!$        print*, "HYPRE HYBRID (AMG)Num Iterations = ", num_iterations 
        
        call HYPRE_ParCSRHybridGetNumIterati (gr_hypreSolver, num_iterations, ierr)
!!$        print*, "HYPRE HYBRID (DSC)Num Iterations = ", num_iterations         

        call parseSolveStatus
        
        call HYPRE_ParCSRHybridGetFinalRelat (gr_hypreSolver, final_res_norm, ierr) 
!!$        print*, "HYPRE HYBRID Relative Residual Norm = ", final_res_norm
        
     end select
     
     if (.NOT. convergedAndOk) then     

        call blabber
        call postBlabber
        
     else if (.NOT. final_res_norm .LE. gr_hypreRelTol) then     

        call blabber
        call postBlabber
        
     else if (.NOT. converged) then     

        print*,"[gr_hypreSolve]: Nonconvergence in subroutine ", & 
             "gr_hypreSolve after",num_iterations," iterations, final_res_norm =",final_res_norm
!!$        print*,' -- converged=',converged,iConverged
        
     else if (gr_hyprePrintSolveInfo) then
        
#if (0)
        print*, "HYPRE SOLVE: Num Iterations = ", num_iterations 
        print*, "HYPRE SOLVE: Relative Residual Norm = ", final_res_norm          
#else
        call blabber
#endif
        hypreError = 0
        call postBlabber
        
     else
        hypreError = 0
        call postBlabber
     end if

  end if


  call HYPRE_ClearAllErrors(ierr)
  
  !! Rob Falgout recommended on 2013-01-09 to do the solver Destroy()
  !! calls explicitly each time after the solver has been called.
  call gr_hypreDestroySolver()

  call Timers_stop("gr_hypreSolve") 
  
  return

contains
  subroutine parseSolveStatus
    convergedAndOk = &
         (converged .AND. (ierr==0) .AND. (num_iterations < gr_hypreMaxIter))
  end subroutine parseSolveStatus

  subroutine blabber
    use Logfile_interface, ONLY: Logfile_stamp
    character(len=32), dimension(7,2) :: strings
    character(len=MAX_STRING_LENGTH) :: hypreErrDescription
    integer, parameter :: d=2
    integer :: i,j,n

    i = 1 
    strings(i,1) = 'ierr'; write(strings(i,2),*) ierr
    call HYPRE_GetError(hypreError)
    if (hypreError .NE. 0) then
       call Hypre_describeerror(hypreError,hypreErrDescription)
       call removeNullChar(hypreErrDescription)
       i = i+1
       strings(i,1) = 'descr'; strings(i,2) = hypreErrDescription
    end if
    i = i+1
    if (component==3) then
       strings(i,1) = 'group no'; write(strings(i,2),*) group
    else
       strings(i,1) = 'component'; write(strings(i,2),*) component
    end if

    if (gr_meshMe .NE. 0) return

    i = i+1
    strings(i,1) = 'converged'; write(strings(i,2),*) converged
    if (iConverged .NE. 0 .AND. iConverged .NE. 1 .AND. &
         iConverged .NE. -2 .AND. gr_hypreSolverType == HYPRE_PCG) then
       i = i+1
       strings(i,1) = 'iConverged'; write(strings(i,2),*) iConverged
    end if
    if (.TRUE. .OR. gr_hypreSolverType == HYPRE_PCG) then
       i = i+1
       strings(i,1) = '|initial RHS|^2'; write(strings(i,2),'(1P,G12.4)') initial_rhs_norm
    end if
    i = i+1
    strings(i,1) = 'num_iterations'; write(strings(i,2),*) num_iterations
    i = i+1
    strings(i,1) = 'final_res_norm'; write(strings(i,2),*) final_res_norm

    n = i
333 format(A, 7(A,'=',A : ', '))
    if (convergedAndOk) then
       print 333,"[gr_hypreSolve]: Ok ", &
            ((trim(adjustl(strings(i,j))),j=1,d), i=1,n)
       call Logfile_stamp(strings(1:n,:), n, d, &
            '[gr_hypreSolve]: Success')
    else
       print 333,"[gr_hypreSolve]: Nonconv./failure ", &
            ((trim(adjustl(strings(i,j))),j=1,d), i=1,n)
       call Logfile_stamp(strings(1:n,:), n, d, &
            '[gr_hypreSolve]: Nonconv./failure')
    end if

  end subroutine blabber

  subroutine postBlabber

    if(component==3) then
       call RadTrans_getDbgContextPtr(dbgContextRt)
       dbgContextRt%component    = component
       dbgContextRt%group        = group
       dbgContextRt%retriable    = 0
       dbgContextRt%libErrCode   = ierr
       if (ierr==0) dbgContextRt%libErrCode = hypreError
       dbgContextRt%flashErrCode = 0
       if (.NOT.convergedAndOk) then
          dbgContextRt%flashErrCode = 1
       else if (dbgContextRt%libErrCode.NE.0) then
          dbgContextRt%flashErrCode = 1
       else if (.NOT. final_res_norm .LE. gr_hypreRelTol) then
          if (final_res_norm .LE. 0.01) then
             dbgContextRt%flashErrCode = 3
          else
             dbgContextRt%flashErrCode = 2
          end if
       end if
       if (dbgContextRt%willingToRetry) then
          if (dbgContextRt%flashErrCode == 1) then
             dbgContextRt%retriable    = 1
          end if
       end if
    end if
  end subroutine postBlabber

end subroutine gr_hypreSolve
