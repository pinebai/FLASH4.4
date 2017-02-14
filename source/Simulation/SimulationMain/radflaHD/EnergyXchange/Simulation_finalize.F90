!!****if* source/Simulation/SimulationMain/radflaHD/EnergyXchange/Simulation_finalize
!!
!! NAME
!!  Simulation_finalize
!!
!! SYNOPSIS
!!
!!  Simulation_finalize()
!!
!! DESCRIPTION
!!
!!  This dummy function cleans up the Simulation unit, deallocates memory, etc.
!!  However, as nothing needs to be done, only this stub is included.
!!
!! ARGUMENTS
!!
!!
!!
!!***

subroutine Simulation_finalize()
  use Simulation_data, ONLY: sim_fileUnitT, sim_fileUnitE

  implicit none

  close(sim_fileUnitT)
  close(sim_fileUnitE)

  return

end subroutine Simulation_finalize
