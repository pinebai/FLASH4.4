!!****if* source/Grid/GridSolvers/Multigrid/gr_hgSolve
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


subroutine gr_hgSolve(gr_iSource, gr_iSoln, gr_iSls, gr_iCorr, SolveBlock, bndTypes, src_fact, dt, chi)

  !=================================================================

  use Grid_data, ONLY: gr_meshMe
  use gr_hgData, ONLY: gr_hgBndTypes, gr_hgSolnIndex, gr_hgMeshRefineMax, &
       gr_hgMaxCorrections, gr_hgPrintNorm, gr_hgMaxResidualNorm, gr_hgAvgSource

  use Timers_interface, ONLY: Timers_start, Timers_stop
  use RuntimeParameters_interface, ONLY: RuntimeParameters_get
  use Logfile_interface, ONLY: Logfile_stamp
  use gr_hgInterface, ONLY: gr_hgBndry, gr_hgInitSource, gr_hgLevelAdd, gr_hgLevelAddScalar, &
       gr_hgLevelMultiplyScalar, gr_hgNorm, gr_hgProlongBndries, gr_hgResidual,                &
       gr_hgRestoreNodetypes, gr_hgRestrict, gr_hgSetZeroBoundary, gr_hgSetMaxLevel,             &
       gr_hgSolveLevel
  use Timers_interface, ONLY: Timers_start, Timers_stop

  implicit none

#include "constants.h"  ! for MASTER_PE
#include "Flash.h"
#include "Multigrid.h"

  integer, intent(in) :: gr_iSource, gr_iSoln, gr_iSls, gr_iCorr
  integer, intent(in) :: bndTypes(6)
  real,intent(IN),OPTIONAL :: dt, chi, src_fact
  external               SolveBlock

  integer             :: n, m, ierr, extrap
  real                :: norm_lhs, norm_rhs, norm_tmp

  character*70        :: internalFile

  !=====================================================================

  !This routine initializes the multigrid and sets up data
  !arrays for data exchange operations for the different levels
  !in the PARAMESH tree hierarchy.  Not needed in cautious mode.

  !Do nothing if we are using PM2 or "cautious" mode.
#ifdef FLASH_GRID_PARAMESH2
#else

#ifdef CAUTIOUS_PARAMESH_CALL
#else
  call Timers_start("amr_morton_process")
  call amr_mg_init()
  call Timers_stop("amr_morton_process")
#endif

#endif


  ! Record current state of some PARAMESH arrays, before we start messing around with them - KW
  call gr_hgRecordNodeTypeState()

  ! Initialize source data.
  gr_hgBndTypes = bndTypes
  gr_hgSolnIndex = gr_iSoln

  call gr_hgInitSource(gr_iSource, gr_iSoln)  !oK
  call gr_hgNorm(0, 2, gr_iSource, norm_rhs, MG_NODES_LEAF_ONLY) !oK

  ! Coarse grid interpolation step.  Solve on coarse grid, interpolate boundary

  ! conditions for next finer level, solve on that level, etc.

  do m = 1, gr_hgMeshRefineMax
     call gr_hgSetMaxLevel(m)    ! This is a new FLASH3 routine
     ! call gr_hgNorm(0, 1, gr_iSource, norm_tmp, MG_NODES_LEAF_ONLY)
     call gr_hgSetZeroBoundary(m, gr_iSoln)  !oK
     call gr_hgSolveLevel(m, gr_iSource, gr_iSoln, SolveBlock, MG_NODES_ALL_NODES)  ! Working here.....
     call gr_hgResidual(m, gr_iSource, gr_iSoln, gr_iSls) !oK
     call gr_hgProlongBndries(m, gr_iSoln, gr_iSoln, 0)  ! LBR doesn't wanna know.  But we looked at it...
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
     do m = gr_hgMeshRefineMax-1, 1, -1
        call gr_hgRestrict(m+1, gr_iSls, gr_iSls)
     enddo
     call Timers_stop("gr_hgRestrict")

     call gr_hgSetMaxLevel(1)
     call gr_hgSetZeroBoundary(1, gr_iCorr)
     call gr_hgSolveLevel(1, gr_iSls, gr_iCorr, SolveBlock, MG_NODES_ALL_NODES, dt, chi)
     call gr_hgLevelAdd(1, gr_iSoln, gr_iCorr, MG_NODES_LEAF_ONLY) !oK
     call gr_hgProlongBndries(1, gr_iCorr, gr_iCorr, 0)  !infamous

     do m = 2, gr_hgMeshRefineMax
        call gr_hgSetMaxLevel(m)
        call gr_hgSetZeroBoundary(m, gr_iCorr)
        call gr_hgSolveLevel(m, gr_iSls, gr_iCorr, SolveBlock, MG_NODES_ALL_NODES)
        call gr_hgResidual(m, gr_iSls, gr_iCorr, gr_iSls)  ! not called on top level
        call gr_hgLevelAdd(m, gr_iSoln, gr_iCorr, MG_NODES_LEAF_ONLY)
        call gr_hgProlongBndries(m, gr_iCorr, gr_iCorr, 0)
     enddo

  enddo

  if (present(src_fact)) then
  ! Multiply solution by the source term factor.
     do m = 1, gr_hgMeshRefineMax
        call gr_hgLevelMultiplyScalar(m, gr_iSoln, src_fact, MG_NODES_LEAF_ONLY) !oK
     enddo
  end if

  ! Add back source average to source function (to compensate for its subtraction
  ! when using periodic/Neumann boundary conditions).

  do m = 1, gr_hgMeshRefineMax
     call gr_hgLevelAddScalar(m, gr_iSource, gr_hgAvgSource, MG_NODES_LEAF_ONLY) !oK
  enddo

  ! Leave boundary zones properly updated.. - SHOULD NOT BE NEEDED! - KW
!!$  call gr_hgBndry (0, gr_iSoln, NGUARD, MG_NODES_LEAF_ONLY, MG_UPDATE_UNK, &
!!$       MG_STANDALONE, .true.)

  ! Make sure we leave PARAMESH in a sane state.
  call gr_hgRestoreNodetypes()
  !==============================================================================

  return
end subroutine gr_hgSolve
