!!****if* source/Grid/GridSolvers/BiPCGStab/poisson/Grid_solvePoisson
!!
!! NAME
!!
!!  Grid_solvePoisson
!!
!! SYNOPSIS
!!
!!  call Grid_solvePoisson(integer(IN) :: iSoln,
!!                         integer(IN) :: iSrc, 
!!                         integer(IN) :: bcTypes(6),
!!                         real(IN)    :: bcValues(2,6),
!!                         real(INOUT) :: poisFact)
!!
!! DESCRIPTION
!!
!!   Driver routine for the BiPCGStab Poisson solver.
!!
!!
!! ARGUMENTS
!!
!!  iSoln    - the index for the solution variable. The solution is
!!             written directly into this variable.
!!  iSrc     - the index of the source variable (density when used for self-gravity)
!!  bcTypes  - the boundary condition type; only the first entry is used.
!!             Only the first 2*NDIM elements are significant. They are interpreted
!!             in the order (X left, X right, Y left, Y right, Z left, Z right).
!!             Valid values are:
!!               GRID_PDE_BND_PERIODIC (1)
!!               GRID_PDE_BND_DIRICHLET (2) (homogeneous or constant Dirichlet)
!!               GRID_PDE_BND_NEUMANN (3) (homogeneous or constant Neumann)
!!               GRID_PDE_BND_ISOLATED (0) - NOT IMPLEMENTED HERE
!!  bcValues - Fixed values at boundary for Dirichlet and Neumann boundary
!!                         conditions.  If Robin, then bcValues(1,*) holds
!!                         Dirichlet and bcValues(2,*) Neumann conditions; if
!!                         Neumann or Dirichlet, then bcValues(1,*) holds
!!                         Neumann or Dirichlet values and bcValues(2,*) goes unused
!!  poisFact      - Constant Poisson factor.  Used to scale the
!!                  source function.  For example, for gravity,
!!                  poisFact = 4*pi*G.
!!
!! NOTES
!!
!!  The symbols listed above for bcTypes are declared as FORTRAN PARAMETERS in
!!  the module Grid_interfaces.  Code using this interface should refer to that
!!  module with a USE statement, like this:
!!
!!    use Grid_interface, ONLY : GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
!!       GRID_PDE_BND_ISOLATED, GRID_PDE_BND_DIRICHLET, &
!!       Grid_solvePoisson
!!  
!!  Most implementations only support a limited subset of boundary condition types.
!!  Some implementations require that all significant elements of bcTypes are the same.
!!  (That is always the case for GRID_PDE_BND_ISOLATED.)
!!  Even if an implementation supports combinations of different boundary conditions
!!  on different sides of the domain, the types at left and right sides for the same
!!  axis direction will usually have to be the samme.
!!
!!  Support in some implementations provided with FLASH4:
!!
!!   GridSolvers/Multipole:                  GRID_PDE_BND_ISOLATED                   
!!   GridSolvers/Multigrid (simple):         GRID_PDE_BND_PERIODIC, GRID_PDE_BND_DIRICHLET(hom.),
!!                                            GRID_PDE_BND_NEUMANN(hom.)(?)
!!                                            (same type in all directions),
!!                                            or GRID_PDE_BND_ISOLATED
!!                                           (requires Paramesh as Grid with NBlockX==NBlockY==NBlockZ==1)
!!   GridSolvers/Pfft:                       GRID_PDE_BND_PERIODIC, GRID_PDE_BND_DIRICHLET(hom.), 
!!                                            GRID_PDE_BND_NEUMANN(hom.)
!!                                            in various combinations,
!!                                            depending on GridSolvers/Pfft subdirectory 
!!                                            (i.e. implementation configured in)
!!                                           (requires UG in pencil shape or Paramesh as Grid)
!!   GridSolvers/Multigrid hybrid with Pfft: GRID_PDE_BND_PERIODIC, GRID_PDE_BND_DIRICHLET, 
!!                                            GRID_PDE_BND_NEUMANN in various combinations,
!!                                            or GRID_PDE_BND_ISOLATED
!!                                           (requires Paramesh as Grid)
!!   GridSolvers/BiPCGStab preconditioned with MC multigrid:
!!                                            GRID_PDE_BND_PERIODIC, GRID_PDE_BND_DIRICHLET(hom.),
!!                                            GRID_PDE_BND_NEUMANN(hom.) in various combinations
!!                                           (requires Paramesh as Grid)
!!   GridSolvers/BiPCGStab preconditioned with (HG) multigrid:
!!                                            GRID_PDE_BND_PERIODIC, GRID_PDE_BND_DIRICHLET(hom.),
!!                                            GRID_PDE_BND_NEUMANN(hom.) in various combinations
!!                                           (requires Paramesh as Grid)
!!   GridSolvers/BiPCGStab preconditioned with HYPRE:
!!                                            Status unknown
!!   
!!***
!*******************************************************************************

!  Routine:     Grid_solvePoisson()

!  Description: Driver routine for the BiPCGStab Poisson solver.  

!  Parameters:  iSrc            Index for source array.  This is taken to be
!                               the density field; the source array to be used
!                               as the right-hand side of the Poisson equation
!                               is computed from this.
!               iSoln           Index for solution array.  The solution is
!                               written directly into this variable.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NOTE: The following have been added to the interface only to get this
! version of the paramesh3.x poisson solver to work.  These have no effect. 
!               bcTypes(6)     Boundary condition types array:
!                                 GRID_PDE_BND_PERIODIC
!                                 GRID_PDE_BND_DIRICHLET
!                                 GRID_PDE_BND_NEUMANN
!
!                                 index 1 = -x, 2 = +x, 3 = -y, 4 = +y, 5 = -z  6 = +z
!               bcValues(2,6)  Values for dirichlet and neumann boundary
!                               conditions.  If Robin, then bcValues(1,*) holds
!                               dirichlet and bcValues(2,*) neumann conditions; if
!                               neumann or dirichlet, then bcValues(1,*) holds
!                               neumann or dirichlet values and bcValues(2,*) goes unused
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!                                 1 = Periodic boundaries
!                                 2 = Dirichlet boundaries
!                                 3 = Neumann boundaries
!               poisfact        Constant Poisson factor.  Used to scale the
!                               source function.  For example, for gravity,
!                               poisfact = 4*pi*G.


subroutine Grid_solvePoisson (iSoln, iSrc, bcTypes, bcValues, poisFact)

!===============================================================================
#include "Flash.h"

use gr_bicgData

use Grid_interface,    ONLY : GRID_PDE_BND_ISOLATED, &
                              GRID_PDE_BND_PERIODIC, &
                              GRID_PDE_BND_DIRICHLET,&
                              GRID_PDE_BND_NEUMANN, &
                              Grid_getLocalNumBlks, &
                              Grid_getListOfBlocks, &
                              Grid_getBlkPtr, &
                              Grid_releaseBlkPtr

!use gr_bicgInterface, ONLY: gr_bipcgstab
use Timers_interface, ONLY : Timers_start, Timers_stop
use Driver_interface, ONLY : Driver_abortFlash


implicit none

  integer, intent(in)    :: iSoln, iSrc
  integer, intent(in)    :: bcTypes(6)
  real, intent(in)       :: bcValues(2,6)
  real, intent(inout)    :: poisFact

integer       :: lb, lnblocks2, MyPE2, MasterPE2
integer       :: i, j, k

integer, save :: bcro,bcri,bcvi,bcpi,bcsi,bczi,bcyi


!===============================================================================

! Get key numbers from the database for the temporary variables we need.

bcro = BIRO_VAR
bcri = BIRI_VAR
bcvi = BIVI_VAR
bcpi = BIPI_VAR
bcsi = BISI_VAR
bczi = BIZI_VAR
bcyi = BIYI_VAR


!===============================================================================

! Call the multigrid Poisson solver for different types of boundary conditions.

#ifdef IMGP_VAR

! This is a hack to make this version of poisson.F90's interface
! consistent with 'blessed' version.


select case (bcTypes(1))

!-------------------------------------------------------------------------------
  print*,'isolated'
  case (GRID_PDE_BND_ISOLATED) ! isolated boundary conditions

    call Driver_abortFlash("[poisson]  Isolated BCs are not supported on BPCGStab.")

!-------------------------------------------------------------------------------

  case (GRID_PDE_BND_PERIODIC) ! periodic boundary conditions

    call gr_bipcgstab (iSrc,iSoln,poisFact, bcro, bcri, bcvi, bcpi, bcsi, bczi, bcyi, &
                    bcTypes,bcValues)

!-------------------------------------------------------------------------------

  case (GRID_PDE_BND_DIRICHLET) ! Dirichlet boundary conditions

    call gr_bipcgstab (iSrc,iSoln,poisFact, bcro, bcri, bcvi, bcpi, bcsi, bczi, bcyi, &
                    bcTypes,bcValues)

!-------------------------------------------------------------------------------

  case (GRID_PDE_BND_NEUMANN) ! Neumann boundary conditions

    call gr_bipcgstab (iSrc,iSoln,poisFact, bcro, bcri, bcvi, bcpi, bcsi, bczi, bcyi, &
                    bcTypes,bcValues)

!-------------------------------------------------------------------------------

  case default
    call Driver_abortFlash("[poisson]  invalid boundary condition type!")


!-------------------------------------------------------------------------------

end select


#else

    call Timers_start("BiPCGStab_solve")

    call gr_bipcgstab (iSrc,iSoln,poisFact, bcro, bcri, bcvi, bcpi, bcsi, bczi, bcyi, &
                    bcTypes,bcValues)

    call Timers_stop("BiPCGStab_solve")
    
#endif


!===============================================================================

return
end subroutine Grid_SolvePoisson
