!!****if* source/Grid/GridSolvers/gr_solversFinalize
!!
!! NAME
!!
!!  gr_solversFinalize
!!
!! 
!! SYNOPSIS
!!
!!  call gr_solversFinalize()
!!
!!
!! DESCRIPTION
!!
!!  This routine finalizes all solvers; ones in 
!!  use will have something more interesting than a stub
!!  compiled in.
!!
!!***

subroutine gr_solversFinalize()

  use Grid_data, ONLY : gr_isolatedBoundaries
  use gr_pfftInterface, ONLY : gr_pfftFinalize
  use gr_hgInterface, ONLY : gr_hgFinalize, gr_hgPfftFinalize
  use gr_mgInterface, ONLY : gr_mgFinalize
  use gr_bicgInterface, ONLY : gr_bicgFinalize
  use gr_bhInterface, ONLY : gr_bhFinalize

  implicit none 

  ! parallel FFT
  call gr_pfftFinalize()

  ! Multipole
  call gr_mpoleFinalize()

  ! Multigrid
  call gr_hgFinalize()
  call gr_mgFinalize()

  ! BiPCGStab
  call gr_bicgFinalize()
 
  ! Isolated multipole
  if (gr_isolatedBoundaries) then
     call gr_isoMpoleFinalize()
  end if

  ! Parallel FFT at coarse level in Multigrid solver.
  call gr_hgPfftFinalize()


  call gr_hypreFinalize()

  ! BHTree
  call gr_bhFinalize()

  return
end subroutine gr_solversFinalize
