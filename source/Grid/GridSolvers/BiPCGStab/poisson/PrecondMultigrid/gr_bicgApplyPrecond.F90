!!****if* source/Grid/GridSolvers/BiPCGStab/poisson/PrecondMultigrid/gr_bicgApplyPrecond
!!
!! NAME
!!  gr_bicgApplyPrecond
!!
!! SYNOPSIS
!!
!!  call gr_bicgApplyPrecond(integer(IN) :: iSoln,
!!                         integer(IN) :: iSrc, 
!!                         integer(IN) :: bcTypes(6),
!!                         real(IN)    :: bcValues(2,6))
!!
!! DESCRIPTION
!!

!!  
!! ARGUMENTS
!!
!!  iSoln    - the index for the solution variable 
!!  iSrc     - the index of the source variable 
!!  bcTypes  - the boundary condition type; only the first entry is used.
!!                    the following conditions are presently supported:
!!                    GRID_PDE_BND_PERIODIC
!!                    GRID_PDE_BND_ISOLATED
!!                    GRID_PDE_BND_DIRICHLET
!!                    GRID_PDE_BND_NEUMANN could most likely be easily implemented
!!  bcValues - an unused argument used to keep the interface standard
!!
!! RESULT
!!  
!!  The variable at iSoln contains the solution to L(u) = rho for the
!!  given boundary value setup
!!
!! PARAMETERS
!!  IMGM_VAR -- grid variable to store the image mass
!!  IMGP_VAR -- grid variable to store the image potential
!!  ISLS_VAR -- residual grid variable
!!  ICOR_VAR -- correction grid variable
!!
!! NOTES
!!
!!  It is currently assumed in this implementation
!!  DEV: But not checked! - KW
!!  that boundary condition types at all sides of the domain are the same. 
!!  What this means to the user  is that only the first value in bcTypes is
!!  checked here, and may be assumed to give the type of boundary condition
!!  for all 2*NDIM directions.
!!  Presently, the FFT package used only supports the same transform done in
!!  all directions.  Using something different, such as PFFT as the top level,
!!  should remedy this.  
!!  The all-sides-the-same boundary values could be for the sake of the
!!  top-level FFT solve (but the dummy argument is not actually used in this
!!  version).
!!
!! SEE ALSO
!! 
!!  gr_hgSolve, gr_isoImageMass
!!  
!!
!!***

subroutine gr_bicgApplyPrecond (iSoln, iSrc, bcTypes, bcValues)

  use gr_hgData, ONLY : gr_hgMeshRefineMax, gr_hgMaxCorrections, &
                        gr_hgMaxResidualNorm, gr_hgPrintNorm,    &
                        gr_hgBndTypes 

  use Grid_data, ONLY : gr_meshMe

!  use gr_isoInterface, ONLY: gr_isoImageMass, gr_isoImageBdry
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_hgInterface, ONLY: gr_hgLevelAdd, gr_hgLevelZero, gr_hgSolve

  use Grid_interface,    ONLY : GRID_PDE_BND_DIRICHLET, &
                                GRID_PDE_BND_NEUMANN,   &  
                                GRID_PDE_BND_PERIODIC
                                
  use gr_bicgData, ONLY : gcellflg

  implicit none

#include "Multigrid.h"
#include "Flash.h"
#include "constants.h"

  integer, intent(IN)               :: iSoln, iSrc
  integer, dimension(6), intent(IN) :: bcTypes
  real, dimension(2,6), intent(IN)  :: bcValues

  real                              :: norm

  integer               :: eachboundary
  real                  :: poisfact 
  logical, save :: first_call = .true.

  integer, dimension(6), save :: bc_types

  external              gr_hgPoissonSolveBlock

  !=======================================================================

! Map GRID_PDE_ boundary conditions to MG_BND_ ones:
  if (first_call) then
  do eachboundary = 1,2*NDIM
   select case(bcTypes(eachBoundary))
     case (GRID_PDE_BND_PERIODIC)
        bc_types(eachBoundary)  = MG_BND_PERIODIC
     case (GRID_PDE_BND_NEUMANN)
        bc_types(eachBoundary)  = MG_BND_NEUMANN
     case (GRID_PDE_BND_DIRICHLET)
        bc_types(eachBoundary)  = MG_BND_DIRICHLET
     case default
           write(*,*) 'In HG Multigrid gr_bicgApplyPrecond: Boundary Condition not supported'
           call Driver_abortFlash('BCs unsupported')
     end select
  enddo

  !gr_hgBndTypes(1:2*NDIM) = bc_types(1:2*NDIM)

  gr_hgMaxCorrections  = -1
  gr_hgMaxResidualNorm = 1.e-5 
  gr_hgPrintNorm = .true.  

  first_call = .false.
  endif

  poisfact = 1.

  call Timers_start("Precond_Multigrid_solve")

  call gr_hgSolve(iSrc, iSoln, ISLS_VAR, ICOR_VAR, gr_hgPoissonSolveBlock, &
                  bc_types, poisfact)


  call Timers_stop("Precond_Multigrid_solve")


  ! No need to fill guardcells before A*x operation in BiPCGStab
  gcellflg = .false.





!!$#ifdef IMGP_VAR
!!$  ! Initializations performed in the gravity unit.
!!$
!!$  ! Call the Huang-Greengard Poisson solver for different types of boundary
!!$  ! conditions.
!!$
!!$  select case (bcTypes(1))
!!$
!!$     !---------------------------------------------------------------
!!$
!!$     ! In order for isolated boundary conditions to work, an appropriate solver implementation
!!$     ! that defines the appropriate mesh variables and implements gr_isoImage{Mass,Bndry}
!!$     ! must be included in the code.
!!$       
!!$  case (MG_BND_ISOLATED)
!!$
!!$     
!!$     call Timers_start("Multigrid_solve")
!!$
!!$     call gr_hgSolve(iSrc, iSoln, ISLS_VAR, ICOR_VAR, gr_hgPoissonSolveBlock, &
!!$          (/(MG_BND_DIRICHLET,i=1,6)/), poisfact)
!!$
!!$     call Timers_stop("Multigrid_solve")
!!$
!!$     do i = 1, gr_hgMeshRefineMax
!!$        call gr_hgLevelZero(i, IMGP_VAR, MG_NODES_ALL_NODES)
!!$     enddo
!!$
!!$     call Timers_start("Isobnd_Mpole_Solve")
!!$
!!$     call gr_isoImageMass(iSoln, IMGM_VAR)
!!$     call gr_isoImageBdry(IMGM_VAR, IMGP_VAR, 1.)
!!$
!!$     call Timers_stop("Isobnd_Mpole_Solve")
!!$
!!$     call Timers_start("Multigrid_solve")
!!$
!!$     call gr_hgSolve(IMGM_VAR, IMGP_VAR, ISLS_VAR, ICOR_VAR, gr_hgPoissonSolveBlock, &
!!$          (/(MG_BND_GIVENVAL,i=1,6)/), -1.)
!!$     
!!$     call Timers_stop ("Multigrid_solve")
!!$
!!$     do i = 1, gr_hgMeshRefineMax
!!$        call gr_hgLevelAdd(i, iSoln, IMGP_VAR, MG_NODES_LEAF_ONLY)
!!$     enddo
!!$     
!!$     !-------------------------------------------------------------------------------
!!$     
!!$  case (MG_BND_PERIODIC)
!!$     
!!$     call Timers_start("Multigrid_solve")
!!$     call gr_hgSolve(iSrc, iSoln, ISLS_VAR, ICOR_VAR, gr_hgPoissonSolveBlock, &
!!$          bcTypes, poisfact)
!!$
!!$     call Timers_stop("Multigrid_solve")
!!$
!!$     !----------------------------------------------------------------
!!$
!!$  case (MG_BND_DIRICHLET)
!!$
!!$     call Timers_start("Multigrid_solve")
!!$
!!$     call gr_hgSolve(iSrc, iSoln, ISLS_VAR, ICOR_VAR, gr_hgPoissonSolveBlock, &
!!$          bcTypes, poisfact)
!!$
!!$     call Timers_stop("Multigrid_solve")
!!$
!!$     !----------------------------------------------------------------------
!!$
!!$  case default
!!$
!!$     call Driver_abortFlash("[poisson]  invalid boundary condition type!")
!!$
!!$     !-------------------------------------------------------------------------------
!!$
!!$  end select
!!$  !===============================================================================
!!$
!!$#else
!!$
!!$     call Timers_start("Multigrid_solve")
!!$
!!$     call gr_hgSolve(iSrc, iSoln, ISLS_VAR, ICOR_VAR, gr_hgPoissonSolveBlock, &
!!$          bcTypes, poisfact)
!!$
!!$     call Timers_stop("Multigrid_solve")
!!$
!!$#endif

 return
end subroutine gr_bicgApplyPrecond
