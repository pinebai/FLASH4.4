!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/RingCubicPath/Simulation_finalize
!!
!! NAME
!!
!!  Simulation_finalize
!!
!! SYNOPSIS
!!
!!  Simulation_finalize ()
!!
!! DESCRIPTION
!!
!!  Ends the simulation smoothly.
!!
!!***

subroutine Simulation_finalize ()

  use Simulation_data,  ONLY : sim_boxRadius, &
                               sim_radialBox

  implicit none
!
!
!     ...Deallocate the box arrays.
!
!
  deallocate (sim_radialBox)
  deallocate (sim_boxRadius)

  return
end subroutine Simulation_finalize
