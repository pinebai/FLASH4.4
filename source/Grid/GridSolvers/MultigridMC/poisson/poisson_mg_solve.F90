!!****if* source/Grid/GridSolvers/MultigridMC/poisson/poisson_mg_solve
!!
!! NAME
!!
!!  poisson_mg_solve
!!
!! SYNOPSIS
!!
!!  call poisson_mg_solve(:: level,
!!                         :: irhs,
!!                         :: ilhs)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   level : 
!!
!!   irhs : 
!!
!!   ilhs : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

!*******************************************************************************

!  Routine:     mg_solve()

!  Description: Coarse-grid solver for the multigrid module.  This version
!               simply relaxes "to convergence" using mg_relax().

!  Parameters:  level       Level to solve on.
!               irhs        Right-hand side (source) of equation.
!               ilhs        Left-hand side (solution) of equation.  Receives
!                           the solution.


subroutine poisson_mg_solve (level, irhs, ilhs )

!===============================================================================

use gr_mgData

use RuntimeParameters_interface, ONLY : RuntimeParameters_get

implicit none

integer            :: level, irhs, ilhs

integer, save      :: mgrid_solve_max_iter
logical, save      :: first_call = .true.


!===============================================================================

if (first_call) then
  call RuntimeParameters_get('mgrid_solve_max_iter', mgrid_solve_max_iter)
  first_call = .false.
end if

call poisson_mg_relax (level, irhs, ilhs, mgrid_solve_max_iter)

!===============================================================================

return
end
