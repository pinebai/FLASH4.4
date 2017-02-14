!!****if* source/Simulation/SimulationMain/unitTest/Roots/x3Polynomials/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!! SYNOPSIS
!!
!!  Simulation_init ()
!!
!! DESCRIPTION
!!
!!  Initializes the parameters for the Roots x3 polynomials unit test.
!!
!!***

subroutine Simulation_init ()

  use  Simulation_data
  use  RuntimeParameters_interface, ONLY: RuntimeParameters_get

  implicit none

  integer :: ut_getFreeFileUnit
!
!
!    ...Get necessary runtime parameters.
!
!
  call RuntimeParameters_get ("basenm",           sim_baseName )
  call RuntimeParameters_get ('sim_printInfo',    sim_printInfo)
!
!
!    ...Prepare the output file.
!
!
  sim_baseName  = adjustl (sim_baseName)     ! to get ready to use 'trim'
  sim_printUnit = ut_getFreeFileUnit ()
  sim_printName = trim (sim_baseName) // "Printout.dat"

  open (unit = sim_printUnit, file = sim_printName)
!
!
!    ...Set the upper bound of relative accuracies, the x3 polynomial coefficients and the
!       exact roots for comparison.
!
!
  call sim_setMaxRelativeAccuracies  ()
  call sim_setx3PolynomialCoeffs     ()
  call sim_setx3PolynomialExactRoots ()
!
!
!    ...Ready!
!
!
  return
end subroutine Simulation_init
