!!****if* source/Simulation/SimulationMain/unitTest/Roots/x2Polynomials/sim_setMaxRelativeAccuracies
!!
!! NAME
!!
!!  sim_setMaxRelativeAccuracies
!! 
!! SYNOPSIS
!!
!!  call sim_setMaxRelativeAccuracies
!!
!! DESCRIPTION
!!
!!  Sets the maximum allowed relative root accuracy for each x2 polynomial. This is a global
!!  accuracy for each x2 polynomial, in the sense that this maximum applies for all the 2 roots
!!  and all the real and imaginary components of the roots.
!!
!! ARGUMENTS 
!!
!!***

subroutine sim_setMaxRelativeAccuracies ()

  use  Driver_interface, ONLY: Driver_abortFlash

  use  Simulation_data,  ONLY: sim_maxRelativeAccuracy,   &
                               sim_numberOfx2Polynomials
  implicit none

  real :: macheps
!
!
!   ...Check array storage size.
!
!
  if (sim_numberOfx2Polynomials > 16) then
      call Driver_abortFlash ('[sim_setMaxRelativeAccuracies] ERROR: Not enough storage for max accuracies!')
  end if
!
!
!   ...Set the maximum relative accuracy for each x2 polynomial.
!
!
  macheps = epsilon (1.0)

  sim_maxRelativeAccuracy (1)  = macheps
  sim_maxRelativeAccuracy (2)  = macheps
  sim_maxRelativeAccuracy (3)  = macheps
  sim_maxRelativeAccuracy (4)  = macheps
  sim_maxRelativeAccuracy (5)  = macheps
  sim_maxRelativeAccuracy (6)  = macheps
  sim_maxRelativeAccuracy (7)  = macheps
  sim_maxRelativeAccuracy (8)  = macheps
  sim_maxRelativeAccuracy (9)  = macheps
  sim_maxRelativeAccuracy (10) = macheps
  sim_maxRelativeAccuracy (11) = macheps
  sim_maxRelativeAccuracy (12) = macheps
  sim_maxRelativeAccuracy (13) = macheps
  sim_maxRelativeAccuracy (14) = macheps
  sim_maxRelativeAccuracy (15) = macheps
  sim_maxRelativeAccuracy (16) = macheps
!
!
!   ...Ready!
!
!
  return
end subroutine sim_setMaxRelativeAccuracies
