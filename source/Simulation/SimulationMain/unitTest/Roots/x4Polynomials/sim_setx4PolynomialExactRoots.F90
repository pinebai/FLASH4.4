!!****if* source/Simulation/SimulationMain/unitTest/Roots/x4Polynomials/sim_setx4PolynomialExactRoots
!!
!! NAME
!!
!!  sim_setx4PolynomialExactRoots
!! 
!! SYNOPSIS
!!
!!  call sim_setx4PolynomialExactRoots
!!
!! DESCRIPTION
!!
!!  Sets the exact x4 polynomial roots that correspond to the x4 polynomial coefficients.
!!
!!              x4 polynomial = x^4 + Ax^3 + Bx^2 + Cx + D
!!
!! ARGUMENTS 
!!
!!***

subroutine sim_setx4PolynomialExactRoots ()

  use  Driver_interface, ONLY: Driver_abortFlash

  use  Simulation_data,  ONLY: sim_numberOfx4Polynomials, &
                               sim_rootsAnalytical

  implicit none

  integer, parameter :: Re = 1
  integer, parameter :: Im = 2
!
!
!   ...Check array storage size.
!
!
  if (sim_numberOfx4Polynomials > 13) then
      call Driver_abortFlash ('[sim_setx4PolynomialExactRoots] ERROR: Not enough storage for exact roots!')
  end if
!
!
!   ...Set the exact roots.
!
!
  sim_rootsAnalytical (1,Re,1) = 1.0e+9
  sim_rootsAnalytical (2,Re,1) = 1.0e+6
  sim_rootsAnalytical (3,Re,1) = 1.0e+3
  sim_rootsAnalytical (4,Re,1) = 1.0
  sim_rootsAnalytical (1,Im,1) = 0.0
  sim_rootsAnalytical (2,Im,1) = 0.0
  sim_rootsAnalytical (3,Im,1) = 0.0
  sim_rootsAnalytical (4,Im,1) = 0.0

  sim_rootsAnalytical (1,Re,2) = 1.003
  sim_rootsAnalytical (2,Re,2) = 1.002
  sim_rootsAnalytical (3,Re,2) = 1.001
  sim_rootsAnalytical (4,Re,2) = 1.000
  sim_rootsAnalytical (1,Im,2) = 0.0
  sim_rootsAnalytical (2,Im,2) = 0.0
  sim_rootsAnalytical (3,Im,2) = 0.0
  sim_rootsAnalytical (4,Im,2) = 0.0

  sim_rootsAnalytical (1,Re,3) = 1.0e+80
  sim_rootsAnalytical (2,Re,3) = -1.0e+74
  sim_rootsAnalytical (3,Re,3) = -1.0e+76
  sim_rootsAnalytical (4,Re,3) = -1.0e+77
  sim_rootsAnalytical (1,Im,3) = 0.0
  sim_rootsAnalytical (2,Im,3) = 0.0
  sim_rootsAnalytical (3,Im,3) = 0.0
  sim_rootsAnalytical (4,Im,3) = 0.0

  sim_rootsAnalytical (1,Re,4) = 1.0e+14
  sim_rootsAnalytical (2,Re,4) = 2.0
  sim_rootsAnalytical (3,Re,4) = 1.0
  sim_rootsAnalytical (4,Re,4) = -1.0
  sim_rootsAnalytical (1,Im,4) = 0.0
  sim_rootsAnalytical (2,Im,4) = 0.0
  sim_rootsAnalytical (3,Im,4) = 0.0
  sim_rootsAnalytical (4,Im,4) = 0.0

  sim_rootsAnalytical (1,Re,5) = 1.0e+7
  sim_rootsAnalytical (2,Re,5) = 1.0
  sim_rootsAnalytical (3,Re,5) = -1.0
  sim_rootsAnalytical (4,Re,5) = -2.0e+7
  sim_rootsAnalytical (1,Im,5) = 0.0
  sim_rootsAnalytical (2,Im,5) = 0.0
  sim_rootsAnalytical (3,Im,5) = 0.0
  sim_rootsAnalytical (4,Im,5) = 0.0

  sim_rootsAnalytical (1,Re,6) = 1.0e+7
  sim_rootsAnalytical (2,Re,6) = -1.0e+6
  sim_rootsAnalytical (3,Re,6) = 1.0
  sim_rootsAnalytical (4,Re,6) = 1.0
  sim_rootsAnalytical (1,Im,6) = 0.0
  sim_rootsAnalytical (2,Im,6) = 0.0
  sim_rootsAnalytical (3,Im,6) = 1.0
  sim_rootsAnalytical (4,Im,6) = -1.0

  sim_rootsAnalytical (1,Re,7) = -4.0
  sim_rootsAnalytical (2,Re,7) = -7.0
  sim_rootsAnalytical (3,Re,7) = -1.0e+6
  sim_rootsAnalytical (4,Re,7) = -1.0e+6
  sim_rootsAnalytical (1,Im,7) = 0.0
  sim_rootsAnalytical (2,Im,7) = 0.0
  sim_rootsAnalytical (3,Im,7) = 1.0e+5
  sim_rootsAnalytical (4,Im,7) = -1.0e+5

  sim_rootsAnalytical (1,Re,8) = 1.0e+8
  sim_rootsAnalytical (2,Re,8) = 11.0
  sim_rootsAnalytical (3,Re,8) = 1.0e+3
  sim_rootsAnalytical (4,Re,8) = 1.0e+3
  sim_rootsAnalytical (1,Im,8) = 0.0
  sim_rootsAnalytical (2,Im,8) = 0.0
  sim_rootsAnalytical (3,Im,8) = 1.0
  sim_rootsAnalytical (4,Im,8) = -1.0

  sim_rootsAnalytical (1,Re,9) = 1.0e+7
  sim_rootsAnalytical (2,Re,9) = 1.0e+7
  sim_rootsAnalytical (3,Re,9) = 1.0
  sim_rootsAnalytical (4,Re,9) = 1.0
  sim_rootsAnalytical (1,Im,9) = 1.0e+6
  sim_rootsAnalytical (2,Im,9) = -1.0e+6
  sim_rootsAnalytical (3,Im,9) = 2.0
  sim_rootsAnalytical (4,Im,9) = -2.0

  sim_rootsAnalytical (1,Re,10) = 1.0e+4
  sim_rootsAnalytical (2,Re,10) = 1.0e+4
  sim_rootsAnalytical (3,Re,10) = -7.0
  sim_rootsAnalytical (4,Re,10) = -7.0
  sim_rootsAnalytical (1,Im,10) = 3.0
  sim_rootsAnalytical (2,Im,10) = -3.0
  sim_rootsAnalytical (3,Im,10) = 1.0e+3
  sim_rootsAnalytical (4,Im,10) = -1.0e+3

  sim_rootsAnalytical (1,Re,11) = 1.002
  sim_rootsAnalytical (2,Re,11) = 1.002
  sim_rootsAnalytical (3,Re,11) = 1.001
  sim_rootsAnalytical (4,Re,11) = 1.001
  sim_rootsAnalytical (1,Im,11) = 4.998
  sim_rootsAnalytical (2,Im,11) = -4.998
  sim_rootsAnalytical (3,Im,11) = 5.001
  sim_rootsAnalytical (4,Im,11) = -5.001

  sim_rootsAnalytical (1,Re,12) = 1.0e+3
  sim_rootsAnalytical (2,Re,12) = 1.0e+3
  sim_rootsAnalytical (3,Re,12) = 1.0e+3
  sim_rootsAnalytical (4,Re,12) = 1.0e+3
  sim_rootsAnalytical (1,Im,12) = 3.0
  sim_rootsAnalytical (2,Im,12) = -3.0
  sim_rootsAnalytical (3,Im,12) = 1.0
  sim_rootsAnalytical (4,Im,12) = -1.0

  sim_rootsAnalytical (1,Re,13) = 2.0
  sim_rootsAnalytical (2,Re,13) = 2.0
  sim_rootsAnalytical (3,Re,13) = 1.0
  sim_rootsAnalytical (4,Re,13) = 1.0
  sim_rootsAnalytical (1,Im,13) = 1.0e+4
  sim_rootsAnalytical (2,Im,13) = -1.0e+4
  sim_rootsAnalytical (3,Im,13) = 1.0e+3
  sim_rootsAnalytical (4,Im,13) = -1.0e+3
!
!
!   ...Ready!
!
!
  return
end subroutine sim_setx4PolynomialExactRoots
