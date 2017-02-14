!!****f* source/Grid/Grid_solvePoisson
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
!!                         real(INOUT) :: poisfact)
!!
!! DESCRIPTION
!!
!!   Driver routine for poisson solvers in the Grid
!!
!!
!! ARGUMENTS
!!
!!  iSoln    - the index for the solution variable (potential when used for self-gravity)
!!  iSrc     - the index of the source variable (density when used for self-gravity)
!!  bcTypes  - the boundary condition type; only the first entry is used.
!!             Only the first 2*NDIM elements are significant. They are interpreted
!!             in the order (X left, X right, Y left, Y right, Z left, Z right).
!!             Valid values are:
!!               GRID_PDE_BND_PERIODIC (1)
!!               GRID_PDE_BND_DIRICHLET (2) (homogeneous or constant Dirichlet)
!!               GRID_PDE_BND_NEUMANN (3) (homogeneous or constant Neumann)
!!               GRID_PDE_BND_ISOLATED (0)
!!  bcValues - the values to boundary conditions, currently not used (treated as 0)
!!  poisfact      - scaling factor to be used in calculation
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
!!   GridSolvers/Multipole_new:              GRID_PDE_BND_ISOLATED                   
!!   GridSolvers/Multigrid (simple):         GRID_PDE_BND_PERIODIC, GRID_PDE_BND_DIRICHLET(hom.),
!!                                            GRID_PDE_BND_NEUMANN(hom.)(?)
!!                                            (same type in all directions),
!!                                            or GRID_PDE_BND_ISOLATED
!!                                           (requires Paramesh as Grid with NBlockX==NBlockY==NBlockZ==1)
!!   GridSolvers/MultigridMC (simple):       GRID_PDE_BND_PERIODIC, GRID_PDE_BND_DIRICHLET(hom.),
!!                                           GRID_PDE_BND_NEUMANN(hom.)
!!                                            (same type in all directions);
!!                                           GRID_PDE_BND_ISOLATED does not work properly.
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
!!   GridSolvers/MultigridMC hybrid with Pfft: GRID_PDE_BND_PERIODIC, GRID_PDE_BND_DIRICHLET, 
!!                                            GRID_PDE_BND_NEUMANN in various combinations;
!!                                            GRID_PDE_BND_ISOLATED does not work properly.
!!                                           (requires Paramesh as Grid)
!!   GridSolvers/BiPCGStab preconditioned with MC multigrid:
!!                                            GRID_PDE_BND_PERIODIC, GRID_PDE_BND_DIRICHLET, 
!!                                            GRID_PDE_BND_NEUMANN in various combinations
!!                                           (requires Paramesh as Grid)
!!   GridSolvers/BiPCGStab preconditioned with (HG) multigrid:
!!                                            GRID_PDE_BND_PERIODIC, GRID_PDE_BND_DIRICHLET, 
!!                                            GRID_PDE_BND_NEUMANN in various combinations
!!                                           (requires Paramesh as Grid)
!!   GridSolvers/BiPCGStab preconditioned with HYPRE:
!!                                            GRID_PDE_BND_PERIODIC, GRID_PDE_BND_DIRICHLET, 
!!                                            GRID_PDE_BND_NEUMANN in various combinations
!!                                           (requires Paramesh as Grid ?)
!!   
!!***

subroutine Grid_solvePoisson (iSoln, iSrc, bcTypes, bcValues, poisfact)


  implicit none

  integer, intent(in)    :: iSoln, iSrc
  integer, intent(in)    :: bcTypes(6)
  real, intent(in)       :: bcValues(2,6)
  real, intent(inout)    :: poisfact !DEV: NOT intent(IN) because some implementation actually changes it? - KW
  
  
  return
end subroutine Grid_solvePoisson
