!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/testII/sim_extractBeamData
!!
!!  NAME 
!!
!!   sim_extractBeamData
!!
!!  SYNOPSIS
!!
!!   sim_extractBeamData ()
!!
!!  DESCRIPTION
!!
!!   This routine extracts some needed beam data for analyzing the outcome of the simulation.
!!
!! ARGUMENTS
!!
!!***

subroutine sim_extractBeamData ()

  use Simulation_data,             ONLY : sim_beamTargetRadius, &
                                          sim_numRaysLaunched,  &
                                          sim_powerPartition

  use ed_extractBeamsData,         ONLY : ed_extractBeamData

  implicit none
!
!
!     ...Retrieve the following data of the beam:
!
!           1) power partition function of the beam grid
!           2) the number of rays that will be launched for the beam
!           3) the beam's radius at the target
!
!
  call ed_extractBeamData   (beamID     = 1,                 &
                             entryField = 'gridWeight',      &
                             dataValue  = sim_powerPartition )

  call ed_extractBeamData   (beamID     = 1,                  &
                             entryField = 'numberOfRays',     &
                             dataValue  = sim_numRaysLaunched )  

  call ed_extractBeamData   (beamID     = 1,                     &
                             entryField = 'targetSemiAxisMajor', &
                             dataValue  = sim_beamTargetRadius   )  
!
!
!     ...Ready!
!
!
  return
end subroutine sim_extractBeamData
