!!****if* source/Grid/GridSolvers/Multigrid/Grid_solvePoisson
!!
!! NAME
!!  Grid_solvePoisson
!!
!! SYNOPSIS
!!
!!  call Grid_solvePoisson(integer(IN) :: iSoln,
!!                         integer(IN) :: iSrc, 
!!                         integer(IN) :: bcTypes(6),
!!                         real(IN)    :: bcValues(2,6),
!!                         real(INOUT) :: poisfact)
!!
!! DESCRIPTION
!!
!!  This routine solves Poisson's for assorted boundary types using the
!!  Huang-Greengard multilevel solver.  It uses the James' image-mass 
!!  method to handle isolated boundaries.
!!  
!! ARGUMENTS
!!
!!  iSoln    - the index for the solution variable (potential when used for self-gravity)
!!  iSrc     - the index of the source variable (density when used for self-gravity)
!!  bcTypes  - the boundary condition type; only the first entry is used to determine
!!             whether all boundaries should be handled as isolated.
!!                    the following conditions are presently supported:
!!                    MG_BND_PERIODIC
!!                    MG_BND_ISOLATED
!!                    MG_BND_DIRICHLET
!!                    MG_BND_NEUMANN could most likely be easily implemented
!!  bcValues - an unused argument used to keep the interface standard
!!  poisfact      - the scaling factor for the eventual solution
!!                  for gravity this is 4*pi*G
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

subroutine Grid_solvePoisson (iSoln, iSrc, bcTypes, bcValues, poisfact)

  use gr_hgData, ONLY : gr_hgMeshRefineMax
  use Grid_data, ONLY : gr_meshMe

  use gr_isoInterface, ONLY: gr_isoImageMass, gr_isoImageBdry
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_hgInterface, ONLY: gr_hgLevelAdd, gr_hgLevelZero, gr_hgSolve

  implicit none

#include "Multigrid.h"
#include "Flash.h"
#include "constants.h"

  integer, intent(IN)               :: iSoln, iSrc
  integer, dimension(6), intent(IN) :: bcTypes
  real, dimension(2,6), intent(IN)  :: bcValues
  real, intent(INOUT)               :: poisfact !DEV: NOT intent(IN) because some implementation actually changes it? - KW
  real                              :: norm

  integer               :: i

  external              gr_hgPoissonSolveBlock

  !=======================================================================


#ifdef IMGP_VAR
  ! Initializations performed in the gravity unit.

  ! Call the Huang-Greengard Poisson solver for different types of boundary
  ! conditions.

  select case (bcTypes(1))

     !---------------------------------------------------------------

     ! In order for isolated boundary conditions to work, an appropriate solver implementation
     ! that defines the appropriate mesh variables and implements gr_isoImage{Mass,Bndry}
     ! must be included in the code.
       
  case (MG_BND_ISOLATED)

     
     call Timers_start("Multigrid_solve")

     call gr_hgSolve(iSrc, iSoln, ISLS_VAR, ICOR_VAR, gr_hgPoissonSolveBlock, &
          (/(MG_BND_DIRICHLET,i=1,6)/), poisfact)

     call Timers_stop("Multigrid_solve")

     do i = 1, gr_hgMeshRefineMax
        call gr_hgLevelZero(i, IMGP_VAR, MG_NODES_ALL_NODES)
     enddo

     call Timers_start("Isobnd_Mpole_Solve")

     call gr_isoImageMass(iSoln, IMGM_VAR)
     call gr_isoImageBdry(IMGM_VAR, IMGP_VAR, 1.)

     call Timers_stop("Isobnd_Mpole_Solve")

     call Timers_start("Multigrid_solve")

     call gr_hgSolve(IMGM_VAR, IMGP_VAR, ISLS_VAR, ICOR_VAR, gr_hgPoissonSolveBlock, &
          (/(MG_BND_GIVENVAL,i=1,6)/), -1.)
     
     call Timers_stop ("Multigrid_solve")

     do i = 1, gr_hgMeshRefineMax
        call gr_hgLevelAdd(i, iSoln, IMGP_VAR, MG_NODES_LEAF_ONLY)
     enddo
     
     !-------------------------------------------------------------------------------
     
  case (MG_BND_PERIODIC)
     
     call Timers_start("Multigrid_solve")
     call gr_hgSolve(iSrc, iSoln, ISLS_VAR, ICOR_VAR, gr_hgPoissonSolveBlock, &
          bcTypes, poisfact)

     call Timers_stop("Multigrid_solve")

     !----------------------------------------------------------------

  case (MG_BND_DIRICHLET)

     call Timers_start("Multigrid_solve")

     call gr_hgSolve(iSrc, iSoln, ISLS_VAR, ICOR_VAR, gr_hgPoissonSolveBlock, &
          bcTypes, poisfact)

     call Timers_stop("Multigrid_solve")

     !----------------------------------------------------------------------

  case default

     call Driver_abortFlash("[poisson]  invalid boundary condition type!")

     !-------------------------------------------------------------------------------

  end select
  !===============================================================================

#else

     call Timers_start("Multigrid_solve")

     call gr_hgSolve(iSrc, iSoln, ISLS_VAR, ICOR_VAR, gr_hgPoissonSolveBlock, &
          bcTypes, poisfact)

     call Timers_stop("Multigrid_solve")

#endif






 return
end subroutine Grid_solvePoisson
