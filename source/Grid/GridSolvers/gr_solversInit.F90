!!****if* source/Grid/GridSolvers/gr_solversInit
!!
!! NAME
!!
!!  gr_solversInit
!!
!! 
!! SYNOPSIS
!!
!!  call gr_solversInit()
!!
!!
!! DESCRIPTION
!!
!!  This routine initializes all grid solvers; solvers in use will
!!  have something more exciting than a stub.
!!
!!***

subroutine gr_solversInit()    

  use Grid_interface, ONLY : Grid_setSolverDbgContextInfo
  use Grid_data, ONLY : gr_isolatedBoundaries
  use gr_pfftInterface, ONLY : gr_pfftInit
  use gr_hgInterface, ONLY : gr_hgInit, gr_hgPfftInit
  use gr_mgInterface, ONLY : gr_mgInit, gr_mgPfftInit
  use gr_bicgInterface, ONLY : gr_bicgInit
  use gr_bhInterface, ONLY : gr_bhInit

  implicit none 

  ! reset debug context info
  call Grid_setSolverDbgContextInfo()

  ! parallel FFT
  call gr_pfftInit()

  ! Multipole
  call gr_mpoleInit()

  ! Multigrid
  call gr_hgInit()
  call gr_mgInit()

  ! BiPCGSTAB
  call gr_bicgInit()

  ! Isolated multipole
  if (gr_isolatedBoundaries) then
     call gr_isoMpoleInit() 
  end if

  ! Parallel FFT at coarse level in Multigrid solver.
  call gr_hgPfftInit()
  call gr_mgPfftInit()

  call gr_hypreInit()

  ! Barnes-Hut Tree
  call gr_bhInit()

  ! solver unit testing
  call gr_solversTestInit()

  return
end subroutine gr_solversInit     
