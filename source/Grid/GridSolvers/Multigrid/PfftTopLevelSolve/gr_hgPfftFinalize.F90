!!****if* source/Grid/GridSolvers/Multigrid/PfftTopLevelSolve/gr_hgPfftFinalize
!!
!! NAME
!!
!!  gr_hgPfftFinalize
!!
!! SYNOPSIS
!!
!!  gr_hgPfftFinalize()
!!
!! DESCRIPTION
!!
!!  This routine cleans up the data necessary for the Multigrid PFFT extensions.
!!
!!***

subroutine gr_hgPfftFinalize()

  use gr_hgInterface, ONLY : gr_hgPfftFinalizeGrid
  implicit none
  
  !When we have completed all the solves we should
  !deallocate any PFFT memory.
  call gr_hgPfftFinalizeGrid()
  
end subroutine gr_hgPfftFinalize
