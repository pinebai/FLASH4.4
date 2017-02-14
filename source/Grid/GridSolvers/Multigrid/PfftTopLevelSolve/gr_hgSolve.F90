!!****if* source/Grid/GridSolvers/Multigrid/PfftTopLevelSolve/gr_hgSolve
!!
!! NAME
!!  gr_hgSolve
!!
!! SYNOPSIS
!!  gr_hgSolve(integer, intent(in) :: gr_iSource,
!!             integer, intent(in) :: gr_iSoln,
!!             integer, intent(in) :: gr_iSls,
!!             integer, intent(in) :: gr_iCorr,
!!             external            :: SolveBlock,
!!             integer, intent(in) :: bndTypes(6),
!!             real(in),OPTIONAL   :: src_fact,
!!             real(in),OPTIONAL   :: dt,
!!             real(in),OPTIONAL   :: chi)
!!
!! DESCRIPTION
!!  This is the main Poisson solve routine for the Huang & Greengard
!!  (2000, SIAM J. Sci. Comput., 21, 1551) algorithm.  This routine 
!!  defines the multigrid cycle as expressed in the article.
!!
!! coarse ^
!!      o---->o                o---->o
!!      ^    s \               ^    s \ c
!!     r|     o o              |r    o o o
!!     e|      l \             |e     l \ r
!!     s|       v o            |s      v o r
!!     t|        e \           |t       e \ e
!!     r|           o          |r          o c
!!     i|            \         |i           \ t
!!     c|             o        |c            o
!!     t|              \       |t             \
!!  src |               o----->o gr_iSls       o------>...until |r| < rtol
!! fine |______________take residual____________take new residual__>
!!
!! ARGUMENTS
!!  gr_iSource - the source variable
!!  gr_iSoln - the solution variable
!!  gr_iSls       - the residual variable
!!  gr_iCorr      - the correction variable
!!  SolveBlock - the name of the function used for individual block solves
!!  bndTypes   - the BC types as defined in Multigrid.h
!!  src_fact   - The factor by which the solution is multiplied
!!  dt         - time step, to be passed down (maybe unused).
!!  chi        - a factor, to be passed down (maybe unused).
!!  
!! EXAMPLE
!!
!!  gr_hgSolve(idens, gpot, hgwk1, hgwk2, gr_hgPoissonSolveBlock,
!!             (/(MG_BND_DIRICHLET,i=1,6)/), newton)
!!
!!  This call to this equation takes the density variable and puts
!!  a 0-boundary-valued solution to poisson's equation in gpot, using
!!  hgwk1 and hgwk2 as work variables (for the residual and correction)
!!  
!!  USED BY
!!
!!   Primarily Grid_solvePoisson.
!!
!!  SEE ALSO
!!   
!!   Grid_solvePoisson
!!   gr_hgSolveLevel
!!    
!!***

#include "constants.h"  ! for MASTER_PE
#include "Flash.h"
#include "Multigrid.h"

subroutine gr_hgSolve(gr_iSource, gr_iSoln, gr_iSls, gr_iCorr, SolveBlock, bndTypes, src_fact, dt, chi)

  use Grid_data, ONLY: gr_meshComm, gr_meshMe, &
       gr_nblockx, gr_nblocky, gr_nblockz
  use gr_hgData, ONLY: gr_hgBndTypes, gr_hgSolnIndex, gr_hgMeshRefineMax, &
       gr_hgMaxCorrections, gr_hgPrintNorm, gr_hgMaxResidualNorm, gr_hgAvgSource
  use gr_hgPfftData, ONLY : gr_hgPfftMaxDirectSolveLevel, gr_hgbcTypes
  use Driver_interface, ONLY : Driver_abortFlash
  use Timers_interface, ONLY: Timers_start, Timers_stop
  use RuntimeParameters_interface, ONLY: RuntimeParameters_get
  use Logfile_interface, ONLY: Logfile_stamp
  use gr_hgInterface, ONLY: gr_hgBndry, gr_hgInitSource, gr_hgLevelAdd, gr_hgLevelAddScalar, &
       gr_hgLevelMultiplyScalar, gr_hgNorm, gr_hgProlongBndries, gr_hgResidual,                &
       gr_hgRestoreNodetypes, gr_hgRestrict, gr_hgSetZeroBoundary, gr_hgSetMaxLevel,             &
       gr_hgSolveLevel, gr_hgPfftInitGrid, gr_hgPfftSolveLevel
  use tree, ONLY : lnblocks, nodetype, grid_changed
  use Grid_interface, ONLY : Grid_getMaxCommonRefinement, &
    Grid_getDeltas
#ifndef FLASH_GRID_PARAMESH2
  use amr_mg_common, ONLY : amr_mg_min_required_level
#endif
  use physicaldata, ONLY : interp_mask_unk
  use workspace, ONLY : interp_mask_work
  implicit none

  integer, intent(in) :: gr_iSource, gr_iSoln, gr_iSls, gr_iCorr
  integer, intent(in) :: bndTypes(6)
  real, optional, intent(IN) :: src_fact, dt, chi
  external               SolveBlock

  integer             :: n, m, eachBoundary
  real                :: norm_lhs, norm_rhs, norm_tmp

  character*70        :: internalFile
  integer             :: comm, globalTopLevel,solveLevel
  logical             :: requestMap
  logical             :: suppressPfft
  integer             :: gridChanged
  integer,save        :: prevSolveLevel = -1

#include "Flash_mpi.h"

  comm = gr_meshComm

  !The PARAMESH variable "grid_changed" takes the value 1 when the 
  !grid has changed.  Otherwise, it has a value 0 when the grid is the same.
  !We store its value locally because if it has value 1, it will be reset to 
  !value 0 in amr_mg_morton_process (called by amr_mg_init). 
  !gr_hgPfftInitGrid needs to know the original value of grid_changed.
  gridChanged = grid_changed
  !=====================================================================
  !write(*,*) 'Interp_mask_unk =',interp_mask_unk
  !write(*,*) 'Interp_mask_work =',interp_mask_work

  !interp_mask_unk = 2  ! 2nd order in Native interpolation gives smaller  
  !interp_mask_work= 2  ! error resp. to analytical sol. (Utest PFFT_POISSON_MULTIGRID)

  ! Determine the maximum common refinement level in the grid.
  call Grid_getMaxCommonRefinement(comm, globalTopLevel)

  !Choose the direct solve level: gr_hgPfftMaxDirectSolveLevel is used to
  !force a solve on a coarser level than globalTopLevel.  We will ignore
  !this value if it is greater than globalTopLevel.
  solveLevel = min(globalTopLevel, gr_hgPfftMaxDirectSolveLevel)

  if (any(bndTypes(1:2*NDIM).GT.MG_BND_GIVENVAL)) then
     solveLevel = 1
     suppressPfft = .TRUE.
     if ((gr_nblockx /= 1) .or. (gr_nblocky /= 1) .or. (gr_nblockz /= 1)) then
        call Driver_abortFlash("[gr_hgSolve] only one block allowed on coarsest level")
     end if
  else
     suppressPfft = .FALSE.
  end if
 
  if (solveLevel < prevSolveLevel) then
     print*,'Lowering solveLevel from',prevSolveLevel,' to',solveLevel,' grid_changed was',grid_changed
     grid_changed = 1
     prevSolveLevel = solveLevel
  end if
  if (prevSolveLevel == -1) prevSolveLevel = solveLevel

  !This routine initializes the multigrid and sets up data
  !arrays for data exchange operations for the different levels
  !in the PARAMESH tree hierarchy.  Not needed in cautious mode.

  !Do nothing if we are using PM2 or "cautious" mode.
#ifdef FLASH_GRID_PARAMESH2
#else

#ifdef CAUTIOUS_PARAMESH_CALL
#else
  call Timers_start("amr_morton_process")
  amr_mg_min_required_level = solveLevel !Set value for amr_mg_init.
  call amr_mg_init()
  call Timers_stop("amr_morton_process")
#endif

#endif

  ! Record current state of some PARAMESH arrays, before we start messing around with them - KW
  call gr_hgRecordNodeTypeState()

  ! Initialize source data.
  ! Initialize PFFT extensions.
  !DEV: CD.  I think we are passing a factor of 1.0 because we are multiplying
  !by the Poisson factor in gr_hgLevelMultiplyScalar, so we don't need to do it twice.
  if (.NOT. suppressPfft) call gr_hgPfftInitGrid(solvelevel, gridChanged, 1.)

  do eachBoundary = 1, 2*NDIM
!!$     if (bndTypes(eachBoundary) == gr_hgbcTypes(eachBoundary) .OR. suppressPfft) then
!!$        gr_hgBndTypes(eachBoundary) = bndTypes(eachBoundary)
     if (bndTypes(eachBoundary) == gr_hgbcTypes(eachBoundary) .OR. suppressPfft .OR. &
          (bndTypes(eachBoundary)==MG_BND_GIVENVAL .AND. gr_hgbcTypes(eachBoundary)==MG_BND_DIRICHLET) ) then
        gr_hgBndTypes(eachBoundary) = bndTypes(eachBoundary)
     else
        if (gr_meshMe .eq. 0) then
           write(*,*) 'gr_hgSolve Error: Boundary Conditions for Poisson Solver is inconsistent.'
           write(*,*) 'gr_hgSolve Error: direction=',eachBoundary
           write(*,*) 'gr_hgSolve Error: gr_hgbcTypes(direction) =',gr_hgbcTypes(eachBoundary)
           write(*,*) 'gr_hgSolve Error: bndTypes(direction)     =',bndTypes(eachBoundary)
        endif
        call Driver_abortFlash('gr_hgSolve Error: BC type argument inconsistent')
     end if
  end do

  gr_hgSolnIndex = gr_iSoln

  call gr_hgInitSource(gr_iSource, gr_iSoln)  !oK
  call gr_hgNorm(0, 2, gr_iSource, norm_rhs, MG_NODES_LEAF_ONLY) !oK

  ! Coarse grid interpolation step.  Solve on coarse grid, interpolate boundary
  ! conditions for next finer level, solve on that level, etc.

  ! --------------------------------------------------------------------------------
  ! We only need to prolong data from solveLevel.  In the original Multigrid 
  ! implementation we had to perform a solve from level 1 to gr_hgMeshRefineMax.  
  ! PFFT saves us the solve from levels 1 to solveLevel-1 because it operates 
  ! on refinement level solveLevel.
  ! --------------------------------------------------------------------------------
  do m = solveLevel, gr_hgMeshRefineMax

     call gr_hgSetMaxLevel(m)                ! This is a new FLASH3 routine
     call gr_hgSetZeroBoundary(m, gr_iSoln)  !oK

     if (m == solveLevel .AND. .NOT. suppressPfft) then
!!$        if (any(bndTypes(1:2*NDIM) == MG_BND_GIVENVAL)) then           
!!$           call fill_outer_sourceCells(solveLevel)
!!$        endif 
        call gr_hgPfftSolveLevel(gr_iSource, gr_iSoln,solveLevel)
     else
        call gr_hgSolveLevel(m, gr_iSource, gr_iSoln, SolveBlock, MG_NODES_ALL_NODES)  ! Working here.....
     end if

     call gr_hgResidual(m, gr_iSource, gr_iSoln, gr_iSls) !oK
     call gr_hgProlongBndries(m, gr_iSoln, gr_iSoln, 0)   ! LBR doesn't wanna know.  But we looked at it...
  enddo

  ! Correction step.  Restrict residuals from finer levels to coarser levels.
  ! Solve for correction on these levels and interpolate boundary conditions to
  ! finer levels.  Solve for corrections there and apply.  Repeat.  Repeat these
  ! correction steps until the desired residual norm is achieved.
  do n = 0, gr_hgMaxCorrections
     call gr_hgNorm(0, 2, gr_iSls, norm_lhs, MG_NODES_LEAF_ONLY)

     if ((gr_hgPrintNorm) .and. (gr_meshMe == MASTER_PE)) then
        write(internalFile, '(A15,I3,A29,ES13.6)') &   ! there was a comma here
             & 'gr_hgSolve: iter ', n, ': norm(residual)/norm(src) = ', norm_lhs/norm_rhs
        write(*,*) internalFile
        call Logfile_stamp(internalFile,'[gr_hgSolve]')
     endif
     if (norm_lhs/norm_rhs <= gr_hgMaxResidualNorm) exit

     call Timers_start("gr_hgRestrict")
     do m = gr_hgMeshRefineMax-1,  solveLevel, -1
        call gr_hgRestrict(m+1, gr_iSls, gr_iSls)
     enddo
     call Timers_stop("gr_hgRestrict")

     call gr_hgSetMaxLevel(solveLevel)
     call gr_hgSetZeroBoundary(solveLevel, gr_iCorr)
     if (.NOT. suppressPfft) then
!!$        if (any(bndTypes(1:2*NDIM) == MG_BND_GIVENVAL)) then           
!!$           call fill_outer_sourceCells(solveLevel)
!!$        endif 
        call gr_hgPfftSolveLevel(gr_iSls, gr_iCorr,solveLevel)
     else
        call gr_hgSolveLevel(1, gr_iSls, gr_iCorr, SolveBlock, MG_NODES_ALL_NODES)
     end if
     if (solveLevel==gr_hgMeshRefineMax) &
          call gr_hgResidual(solveLevel, gr_iSls, gr_iCorr, gr_iSls)  ! not called on top level of several
     call gr_hgLevelAdd(solveLevel, gr_iSoln, gr_iCorr, MG_NODES_LEAF_ONLY) !oK
     call gr_hgProlongBndries(solveLevel, gr_iCorr, gr_iCorr, 0)            !infamous


     ! *** Again we ignore the levels below solveLevel ***
     do m = solveLevel+1, gr_hgMeshRefineMax
        call gr_hgSetMaxLevel(m)
        call gr_hgSetZeroBoundary(m, gr_iCorr)
        call gr_hgSolveLevel(m, gr_iSls, gr_iCorr, SolveBlock, MG_NODES_ALL_NODES)
        call gr_hgResidual(m, gr_iSls, gr_iCorr, gr_iSls)  ! not called on top level of several
        call gr_hgLevelAdd(m, gr_iSoln, gr_iCorr, MG_NODES_LEAF_ONLY)
        call gr_hgProlongBndries(m, gr_iCorr, gr_iCorr, 0)
     enddo

  enddo

  ! Multiply solution by the source term factor.
  if (present(src_fact)) then
    do m = solveLevel, gr_hgMeshRefineMax
      call gr_hgLevelMultiplyScalar(m, gr_iSoln, src_fact, MG_NODES_LEAF_ONLY) !oK
    end do
  end if

  ! Add back source average to source function (to compensate for its subtraction
  ! when using periodic/Neumann boundary conditions).
  do m = solveLevel, gr_hgMeshRefineMax
     call gr_hgLevelAddScalar(m, gr_iSource, gr_hgAvgSource, MG_NODES_LEAF_ONLY) !oK
  enddo

  ! Leave boundary zones properly updated. - SHOULD NOT BE NEEDED! - KW
!!$  call gr_hgBndry (0, gr_iSoln, NGUARD, MG_NODES_LEAF_ONLY, MG_UPDATE_UNK, &
!!$       MG_STANDALONE, .true.)

  ! Make sure we leave PARAMESH in a sane state.
  call gr_hgRestoreNodetypes()
  !==============================================================================

  !interp_mask_unk = 1
  !interp_mask_work= 1

  return

contains
  subroutine fill_outer_sourceCells(level)
    use physicaldata, ONLY : unk
    use tree, ONLY : lrefine
    integer, intent(IN) :: level
    real, dimension(MDIM)        :: deltas
    real, parameter              :: coeff2 = 2.0 !13./6.
    integer :: b
    logical :: SolveThisBlock
    real, pointer ::soln(:,:,:)

    do b = 1, lnblocks

       SolveThisBlock = (lrefine(b) == level)
       SolveThisBlock = (SolveThisBlock .and. (nodetype(b) == LEAF))

       if (SolveThisBlock) then

          !Let's try a standard way of finding delta x rather than this convolution
          call Grid_getDeltas(b,deltas)


          unk(gr_iSource,NGUARD+1,1+NGUARD*K2D:NYB+NGUARD*K2D,1+NGUARD*K3D:NZB+NGUARD*K3D,b) &
                  = unk(gr_iSource,NGUARD+1,1+NGUARD*K2D:NYB+NGUARD*K2D,1+NGUARD*K3D:NZB+NGUARD*K3D,b) &
                  - coeff2*unk(gr_iSoln,NGUARD, &
               1+NGUARD*K2D:NYB+NGUARD*K2D, &
               1+NGUARD*K3D:NZB+NGUARD*K3D,b)/deltas(IAXIS)**2
          unk(gr_iSource,NGUARD+NXB,1+NGUARD*K2D:NYB+NGUARD*K2D,1+NGUARD*K3D:NZB+NGUARD*K3D,b) &
                = unk(gr_iSource,NGUARD+NXB,1+NGUARD*K2D:NYB+NGUARD*K2D,1+NGUARD*K3D:NZB+NGUARD*K3D,b) &
                - coeff2*unk(gr_iSoln,NGUARD+NXB+1, &
               1+NGUARD*K2D:NYB+NGUARD*K2D, &
               1+NGUARD*K3D:NZB+NGUARD*K3D,b)/deltas(IAXIS)**2

          if (NDIM >= 2) then
             unk(gr_iSource,NGUARD+1:NGUARD+NXB,NGUARD+1,1+NGUARD*K3D:NZB+NGUARD*K3D,b) &
                  = unk(gr_iSource,NGUARD+1:NGUARD+NXB,NGUARD+1,1+NGUARD*K3D:NZB+NGUARD*K3D,b) &
                  - coeff2*unk(gr_iSoln,NGUARD+1:NGUARD+NXB, &
                  NGUARD, &
                  1+NGUARD*K3D:NZB+NGUARD*K3D,b)/deltas(JAXIS)**2
             unk(gr_iSource,NGUARD+1:NGUARD+NXB,NGUARD+NYB,1+NGUARD*K3D:NZB+NGUARD*K3D,b) &
                = unk(gr_iSource,NGUARD+1:NGUARD+NXB,NGUARD+NYB,1+NGUARD*K3D:NZB+NGUARD*K3D,b) &
                - coeff2*unk(gr_iSoln,NGUARD+1:NGUARD+NXB, &
                  NGUARD+NYB+1, &
                  1+NGUARD*K3D:NZB+NGUARD*K3D,b)/deltas(JAXIS)**2
          endif

          if (NDIM == 3) then
             unk(gr_iSource,NGUARD+1:NGUARD+NXB,NGUARD+1:NGUARD+NYB,NGUARD+1,b) &
                  = unk(gr_iSource,NGUARD+1:NGUARD+NXB,NGUARD+1:NGUARD+NYB,NGUARD+1,b) &
                  - coeff2*unk(gr_iSoln,NGUARD+1:NGUARD+NXB, &
                  1+NGUARD*K2D:NYB+NGUARD*K2D, &
                  NGUARD,b )/deltas(KAXIS)**2
             unk(gr_iSource,NGUARD+1:NGUARD+NXB,NGUARD+1:NGUARD+NYB,NGUARD+NZB,b) &
                = unk(gr_iSource,NGUARD+1:NGUARD+NXB,NGUARD+1:NGUARD+NYB,NGUARD+NZB,b) &
                - coeff2*unk(gr_iSoln,NGUARD+1:NGUARD+NXB, &
                  1+NGUARD*K3D:NYB+NGUARD*K3D, &
                  (NGUARD+NZB)*K3D+1,b           )/deltas(KAXIS)**2
          endif
       end if
    end do

  end subroutine fill_outer_sourceCells
end subroutine gr_hgSolve
