!!****if* source/Grid/GridSolvers/MultigridMC/poisson/PfftTopLevelSolve/gr_mgPfftFinalize
!!
!! NAME
!!
!!  gr_mgPfftFinalize
!!
!! SYNOPSIS
!!
!!  gr_mgPfftFinalize()
!!
!! DESCRIPTION
!!
!!  This routine cleans up the data necessary for the Multigrid PFFT extensions.
!!
!!***

subroutine gr_mgPfftFinalize()

  use gr_mgInterface, ONLY : gr_mgPfftFinalizeGrid
  implicit none
  
  !When we have completed all the solves we should
  !deallocate any PFFT memory.
  call gr_mgPfftFinalizeGrid()
  
end subroutine gr_mgPfftFinalize
