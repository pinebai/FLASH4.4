!!****if* source/Simulation/SimulationMain/unitTest/Roots/x2Polynomials/sim_setx2PolynomialExactRoots
!!
!! NAME
!!
!!  sim_setx2PolynomialExactRoots
!! 
!! SYNOPSIS
!!
!!  call sim_setx2PolynomialExactRoots
!!
!! DESCRIPTION
!!
!!  Sets the exact x2 polynomial roots that correspond to the x2 polynomial coefficients.
!!
!!                 x2 polynomial = x^2 + Ax + B
!!
!! ARGUMENTS 
!!
!!***

subroutine sim_setx2PolynomialExactRoots ()

  use  Driver_interface, ONLY: Driver_abortFlash

  use  Simulation_data,  ONLY: sim_LPN,                   &
                               sim_numberOfx2Polynomials, &
                               sim_rootsAnalytical,       &
                               sim_sqrtLPN,               &
                               sim_sqrtSPN

  implicit none

  integer, parameter :: Re = 1
  integer, parameter :: Im = 2
!
!
!   ...Check array storage size.
!
!
  if (sim_numberOfx2Polynomials > 16) then
      call Driver_abortFlash ('[sim_setx2PolynomialExactRoots] ERROR: Not enough storage for exact roots!')
  end if
!
!
!   ...Set the exact roots.
!      LL polynomials first.
!
!
  sim_rootsAnalytical (1,Re,1) = - 1.0
  sim_rootsAnalytical (2,Re,1) = - sim_LPN
  sim_rootsAnalytical (1,Im,1) = 0.0
  sim_rootsAnalytical (2,Im,1) = 0.0

  sim_rootsAnalytical (1,Re,2) = 1.0
  sim_rootsAnalytical (2,Re,2) = - sim_LPN
  sim_rootsAnalytical (1,Im,2) = 0.0
  sim_rootsAnalytical (2,Im,2) = 0.0

  sim_rootsAnalytical (1,Re,3) = sim_LPN
  sim_rootsAnalytical (2,Re,3) = 1.0
  sim_rootsAnalytical (1,Im,3) = 0.0
  sim_rootsAnalytical (2,Im,3) = 0.0

  sim_rootsAnalytical (1,Re,4) = sim_LPN
  sim_rootsAnalytical (2,Re,4) = - 1.0
  sim_rootsAnalytical (1,Im,4) = 0.0
  sim_rootsAnalytical (2,Im,4) = 0.0
!
!
!   ...LS polynomials.
!
!
  sim_rootsAnalytical (1,Re,5) = 0.0
  sim_rootsAnalytical (2,Re,5) = - sim_LPN
  sim_rootsAnalytical (1,Im,5) = 0.0
  sim_rootsAnalytical (2,Im,5) = 0.0

  sim_rootsAnalytical (1,Re,6) = 0.0
  sim_rootsAnalytical (2,Re,6) = - sim_LPN
  sim_rootsAnalytical (1,Im,6) = 0.0
  sim_rootsAnalytical (2,Im,6) = 0.0

  sim_rootsAnalytical (1,Re,7) = sim_LPN
  sim_rootsAnalytical (2,Re,7) = 0.0
  sim_rootsAnalytical (1,Im,7) = 0.0
  sim_rootsAnalytical (2,Im,7) = 0.0

  sim_rootsAnalytical (1,Re,8) = sim_LPN
  sim_rootsAnalytical (2,Re,8) = 0.0
  sim_rootsAnalytical (1,Im,8) = 0.0
  sim_rootsAnalytical (2,Im,8) = 0.0
!
!
!   ...SL polynomials.
!
!
  sim_rootsAnalytical (1,Re,9) = 0.0
  sim_rootsAnalytical (2,Re,9) = 0.0
  sim_rootsAnalytical (1,Im,9) = sim_sqrtLPN
  sim_rootsAnalytical (2,Im,9) = - sim_sqrtLPN

  sim_rootsAnalytical (1,Re,10) = sim_sqrtLPN
  sim_rootsAnalytical (2,Re,10) = - sim_sqrtLPN
  sim_rootsAnalytical (1,Im,10) = 0.0
  sim_rootsAnalytical (2,Im,10) = 0.0

  sim_rootsAnalytical (1,Re,11) = 0.0
  sim_rootsAnalytical (2,Re,11) = 0.0
  sim_rootsAnalytical (1,Im,11) = sim_sqrtLPN
  sim_rootsAnalytical (2,Im,11) = - sim_sqrtLPN

  sim_rootsAnalytical (1,Re,12) = sim_sqrtLPN
  sim_rootsAnalytical (2,Re,12) = - sim_sqrtLPN
  sim_rootsAnalytical (1,Im,12) = 0.0
  sim_rootsAnalytical (2,Im,12) = 0.0
!
!
!   ...SS polynomials.
!
!
  sim_rootsAnalytical (1,Re,13) = 0.0
  sim_rootsAnalytical (2,Re,13) = 0.0
  sim_rootsAnalytical (1,Im,13) = sim_sqrtSPN
  sim_rootsAnalytical (2,Im,13) = - sim_sqrtSPN

  sim_rootsAnalytical (1,Re,14) = sim_sqrtSPN
  sim_rootsAnalytical (2,Re,14) = - sim_sqrtSPN
  sim_rootsAnalytical (1,Im,14) = 0.0
  sim_rootsAnalytical (2,Im,14) = 0.0

  sim_rootsAnalytical (1,Re,15) = 0.0
  sim_rootsAnalytical (2,Re,15) = 0.0
  sim_rootsAnalytical (1,Im,15) = sim_sqrtSPN
  sim_rootsAnalytical (2,Im,15) = - sim_sqrtSPN

  sim_rootsAnalytical (1,Re,16) = sim_sqrtSPN
  sim_rootsAnalytical (2,Re,16) = - sim_sqrtSPN
  sim_rootsAnalytical (1,Im,16) = 0.0
  sim_rootsAnalytical (2,Im,16) = 0.0
!
!
!   ...Ready!
!
!
  return
end subroutine sim_setx2PolynomialExactRoots
