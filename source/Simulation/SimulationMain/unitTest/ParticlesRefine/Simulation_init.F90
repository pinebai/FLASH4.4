!!****if* source/Simulation/SimulationMain/unitTest/ParticlesRefine/Simulation_init
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
  
  use Simulation_data, ONLY :  sim_ptMass, sim_densityThreshold,&
       sim_smlRho,  sim_meshMe, sim_minBlks
  use Driver_interface, ONLY : Driver_getMype

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp
  implicit none
#include "Flash.h"
#include "constants.h"

  call Driver_getMype(MESH_COMM, sim_meshMe)

  call Logfile_stamp( "initializing for unit test for particles based refinement",  &
       "[Simulation_init]")
  call RuntimeParameters_get("sim_ptMass",sim_ptMass)
  call RuntimeParameters_get("sim_densityThreshold",sim_densityThreshold)
  call RuntimeParameters_get("sim_smlRho",sim_smlRho)
  call RuntimeParameters_get("sim_minBlks",sim_minBlks)

end subroutine Simulation_init
