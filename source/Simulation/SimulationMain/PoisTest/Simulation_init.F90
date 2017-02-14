!!****if* source/Simulation/SimulationMain/PoisTest/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!!
!! DESCRIPTION
!!
!!  Initializes all the parameters needed for the Huang & Greengard (poistest)
!!  problem
!!
!! ARGUMENTS
!!
!!  
!!
!! PARAMETERS
!!
!!  sim_smlRho : 
!!  sim_gam : 
!!***

subroutine Simulation_init()
  
  use Simulation_data
  use Driver_interface, ONLY : Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp
  implicit none
#include "Flash.h"
#include "constants.h"

  
  call Driver_getMype(MESH_COMM, sim_meshMe)

  call RuntimeParameters_get('sim_smlRho', sim_smlRho)
  
  call Logfile_stamp( "initializing for Huang & Greengard problem",  &
       "[Simulation_init]")

  sim_gam = 1.5

end subroutine Simulation_init
