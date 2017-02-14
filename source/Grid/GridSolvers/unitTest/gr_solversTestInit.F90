!!****if* source/Grid/GridSolvers/unitTest/gr_solversTestInit
!!
!! NAME
!!
!!  gr_solversInit
!!
!! 
!! SYNOPSIS
!!
!!  call gr_solversTestInit()
!!
!!
!! DESCRIPTION
!!
!!  This routine initializes data used by tests
!!  of grid solvers.
!!
!!***

subroutine gr_solversTestInit()
  use RuntimeParameters_interface, ONLY: RuntimeParameters_get
  use gr_solversTestData, ONLY: gr_testTolL2, gr_testTolLinf
  implicit none 

  call RuntimeParameters_get ("gr_testTolL2",   gr_testTolL2)
  call RuntimeParameters_get ("gr_testTolLinf", gr_testTolLinf)

end subroutine gr_solversTestInit     
