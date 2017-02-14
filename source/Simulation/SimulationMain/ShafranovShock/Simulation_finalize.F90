!!****if* source/Simulation/SimulationMain/ShafranovShock/Simulation_finalize
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

  use Simulation_data

  implicit none

  deallocate (sim_x)
  deallocate (sim_velx)
  deallocate (sim_tele)
  deallocate (sim_rho)
  deallocate (sim_tion)
  deallocate (sim_pele)
  deallocate (sim_pion)
  deallocate (sim_sele)

  return

end subroutine Simulation_finalize
