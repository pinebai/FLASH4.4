!!****if* source/Simulation/SimulationMain/FlameChannel/Simulation_sendOutputData
!!
!! NAME
!!  Simulation_sendOutputData
!!
!! SYNOPSIS
!! 
!!  Simulation_sendOutputData()
!!  
!! DESCRIPTION 
!!
!! This routine sends the scalar variables owned by the Simulation unit
!! to the IO unit, to be written to a checkpoint file.
!!
!!
!!***

subroutine Simulation_sendOutputData()

  use Simulation_data, ONLY : sim_last_burned_mass, sim_inflowVx
  use IO_interface, ONLY: IO_setScalar
  
  implicit none

  call IO_setScalar("inflowVx", sim_inflowVx)
  call IO_setScalar("last_burned_mass", sim_last_burned_mass)

end subroutine Simulation_sendOutputData

