!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/RingExpPotential/Simulation_finalize
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
                               sim_boxXvalue, &
                               sim_radialBox, &
                               sim_xvalueBox

  implicit none
!
!
!     ...Deallocate the box arrays.
!
!
  deallocate (sim_radialBox)
  deallocate (sim_xvalueBox)
  deallocate (sim_boxRadius)
  deallocate (sim_boxXvalue)

  return
end subroutine Simulation_finalize
