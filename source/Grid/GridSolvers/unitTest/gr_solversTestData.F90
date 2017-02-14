!!****if* source/Grid/GridSolvers/unitTest/gr_solversTestData
!!
!! NAME
!!
!!  gr_solversTestData
!!
!! SYNOPSIS
!!  use gr_solversTestData
!!
!! DESCRIPTION
!!
!!  Defines some data items that are private to
!!  the GridSolver unit test implementation.
!!
!!  
!!***

Module gr_solversTestData 
  
  implicit none

  ! store error norm tolerances
  real, save :: gr_testTolL2, gr_testTolLinf

end Module gr_solversTestData
