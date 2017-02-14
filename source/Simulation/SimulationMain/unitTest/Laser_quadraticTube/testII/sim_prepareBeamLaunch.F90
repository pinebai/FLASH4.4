!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/testII/sim_prepareBeamLaunch
!!
!!  NAME 
!!
!!   sim_prepareBeamLaunch
!!
!!  SYNOPSIS
!!
!!   sim_prepareBeamLaunch ()
!!
!!  DESCRIPTION
!!
!!   This routine prepares for the launch of 1 circular beam of rays onto the xz-plane of a
!!   quadratic potential tube. It sets the inverse gaussian radial decay radius of the beam equal
!!   to:
!!
!!                   Rx = Ry = sqrt [ (nw + nw) / (A * nuw * tcross) ]
!!
!!   by overriding the value given in the flash.par file. This choice of the inverse gaussian radial
!!   decay radius ensures equal exit power for all rays.
!!
!! ARGUMENTS
!!
!!***

subroutine sim_prepareBeamLaunch ()

  use Simulation_data,             ONLY : sim_A,                &
                                          sim_nuw,              &
                                          sim_nw,               &
                                          sim_tcross

  use ed_overrideBeamsData,        ONLY : ed_overrideBeamData
  use ed_interface,                ONLY : ed_beamsInfo

  implicit none

  real :: Rx, Ry
!
!
!     ...Calculate the gaussian decay radius and store it into the Energy Deposition unit.
!
!
  Rx = sqrt ( (sim_nw + sim_nw) / (sim_A * sim_nuw * sim_tcross) )
  Ry = Rx

  call ed_overrideBeamData  (beamID     = 1,                      &
                             entryField = 'gaussianRadiusMajor', &
                             dataValue  = Rx                      )
  
  call ed_overrideBeamData  (beamID     = 1,                      &
                             entryField = 'gaussianRadiusMinor', &
                             dataValue  = Ry                      )
!
!
!     ...For non-statistical grids we have to recalculate the power partition function, as
!        this property of the beam is permanent and does not change. The statistical grids are
!        treated differently in the sense that the power partition function is calculated on
!        the fly when creating the beam's rays (i.e. when launching the beam).
!
!
  call ed_beamsInfo ()
!
!
!     ...Ready!
!
!
  return
end subroutine sim_prepareBeamLaunch
