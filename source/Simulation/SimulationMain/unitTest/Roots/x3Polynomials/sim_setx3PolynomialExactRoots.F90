!!****if* source/Simulation/SimulationMain/unitTest/Roots/x3Polynomials/sim_setx3PolynomialExactRoots
!!
!! NAME
!!
!!  sim_setx3PolynomialExactRoots
!! 
!! SYNOPSIS
!!
!!  call sim_setx3PolynomialExactRoots
!!
!! DESCRIPTION
!!
!!  Sets the exact x3 polynomial roots that correspond to the x3 polynomial coefficients.
!!
!!              x3 polynomial = x^3 + Ax^2 + Bx^1 + C
!!
!! ARGUMENTS 
!!
!!***

subroutine sim_setx3PolynomialExactRoots ()

  use  Driver_interface, ONLY: Driver_abortFlash

  use  Simulation_data,  ONLY: sim_numberOfx3Polynomials, &
                               sim_rootsAnalytical

  implicit none

  integer, parameter :: Re = 1
  integer, parameter :: Im = 2
!
!
!   ...Check array storage size.
!
!
  if (sim_numberOfx3Polynomials > 9) then
      call Driver_abortFlash ('[sim_setx3PolynomialExactRoots] ERROR: Not enough storage for exact roots!')
  end if
!
!
!   ...Set the exact roots.
!
!
  sim_rootsAnalytical (1,Re,1) = 1.0e+14
  sim_rootsAnalytical (2,Re,1) = 1.0e+7
  sim_rootsAnalytical (3,Re,1) = 1.0
  sim_rootsAnalytical (1,Im,1) = 0.0
  sim_rootsAnalytical (2,Im,1) = 0.0
  sim_rootsAnalytical (3,Im,1) = 0.0

  sim_rootsAnalytical (1,Re,2) = 1.000002
  sim_rootsAnalytical (2,Re,2) = 1.000001
  sim_rootsAnalytical (3,Re,2) = 1.000000
  sim_rootsAnalytical (1,Im,2) = 0.0
  sim_rootsAnalytical (2,Im,2) = 0.0
  sim_rootsAnalytical (3,Im,2) = 0.0

  sim_rootsAnalytical (1,Re,3) = 1.0e+80
  sim_rootsAnalytical (2,Re,3) = 1.0e+77
  sim_rootsAnalytical (3,Re,3) = -1.0e+81
  sim_rootsAnalytical (1,Im,3) = 0.0
  sim_rootsAnalytical (2,Im,3) = 0.0
  sim_rootsAnalytical (3,Im,3) = 0.0

  sim_rootsAnalytical (1,Re,4) = 1.0
  sim_rootsAnalytical (2,Re,4) = -1.0
  sim_rootsAnalytical (3,Re,4) = -1.0e+24
  sim_rootsAnalytical (1,Im,4) = 0.0
  sim_rootsAnalytical (2,Im,4) = 0.0
  sim_rootsAnalytical (3,Im,4) = 0.0

  sim_rootsAnalytical (1,Re,5) = 1.0e+14
  sim_rootsAnalytical (2,Re,5) = 1.0e+14
  sim_rootsAnalytical (3,Re,5) = -1.0
  sim_rootsAnalytical (1,Im,5) = 0.0
  sim_rootsAnalytical (2,Im,5) = 0.0
  sim_rootsAnalytical (3,Im,5) = 0.0

  sim_rootsAnalytical (1,Re,6) = 1.0e+5
  sim_rootsAnalytical (2,Re,6) = 1.0e+5
  sim_rootsAnalytical (3,Re,6) = 1.0e+5
  sim_rootsAnalytical (1,Im,6) = 0.0
  sim_rootsAnalytical (2,Im,6) = 1.0
  sim_rootsAnalytical (3,Im,6) = -1.0

  sim_rootsAnalytical (1,Re,7) = 1.0
  sim_rootsAnalytical (2,Re,7) = 1.0
  sim_rootsAnalytical (3,Re,7) = 1.0
  sim_rootsAnalytical (1,Im,7) = 0.0
  sim_rootsAnalytical (2,Im,7) = 1.0e+7
  sim_rootsAnalytical (3,Im,7) = -1.0e+7

  sim_rootsAnalytical (1,Re,8) = 1.0
  sim_rootsAnalytical (2,Re,8) = 1.0e+7
  sim_rootsAnalytical (3,Re,8) = 1.0e+7
  sim_rootsAnalytical (1,Im,8) = 0.0
  sim_rootsAnalytical (2,Im,8) = 1.0
  sim_rootsAnalytical (3,Im,8) = -1.0

  sim_rootsAnalytical (1,Re,9) = -1.0e+14
  sim_rootsAnalytical (2,Re,9) = 1.0
  sim_rootsAnalytical (3,Re,9) = 1.0
  sim_rootsAnalytical (1,Im,9) = 0.0
  sim_rootsAnalytical (2,Im,9) = 1.0
  sim_rootsAnalytical (3,Im,9) = -1.0
!
!
!   ...Ready!
!
!
  return
end subroutine sim_setx3PolynomialExactRoots
