!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/testII/sim_setErrorBars
!!
!!  NAME 
!!
!!   sim_setErrorBars
!!
!!  SYNOPSIS
!!
!!   sim_setErrorBars ()
!!
!!  DESCRIPTION
!!
!!   This routine sets the error bars for the simulation. These were determined
!!   from previous runs.
!!
!!***

subroutine sim_setErrorBars ()

  use Simulation_data,    ONLY : sim_maxDeltaRsqrErrorBar, &
                                 sim_maxDeltaPZErrorBar

  implicit none
!
!
!     ...Set the upper bounds of the errors. Since these depend on the refinement level,
!        each level has its own upper bound. They also depend on the number of rays launched,
!        so these are for the current specified # of rays in the flash.par files (1,000,000 rays).
!
!
  sim_maxDeltaRsqrErrorBar (1) = 3.4
  sim_maxDeltaRsqrErrorBar (2) = 1.4
  sim_maxDeltaRsqrErrorBar (3) = 3.8e-1
  sim_maxDeltaRsqrErrorBar (4) = 1.3e-1
  sim_maxDeltaRsqrErrorBar (5) = 3.2e-2
  sim_maxDeltaRsqrErrorBar (6) = 8.0e-3

  sim_maxDeltaPZErrorBar   (1) = 1.1e-2
  sim_maxDeltaPZErrorBar   (2) = 9.6e-3
  sim_maxDeltaPZErrorBar   (3) = 5.8e-3
  sim_maxDeltaPZErrorBar   (4) = 4.5e-3
  sim_maxDeltaPZErrorBar   (5) = 2.4e-3
  sim_maxDeltaPZErrorBar   (6) = 1.2e-3
!
!
!     ...Ready!
!
!
  return
end subroutine sim_setErrorBars
