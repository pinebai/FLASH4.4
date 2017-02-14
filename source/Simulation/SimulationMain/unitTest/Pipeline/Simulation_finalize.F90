!!****if* source/Simulation/SimulationMain/unitTest/Pipeline/Simulation_finalize
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

  use Simulation_data, ONLY: sim_items

  implicit none
!
!
!       ...Deallocate the items array.
!
!
  if (allocated (sim_items)) deallocate (sim_items)
!
!
!     ...Ready!
!
!
  return
end subroutine Simulation_finalize
