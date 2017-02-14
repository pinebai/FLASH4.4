!!****if* source/Simulation/SimulationMain/StirTurb/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init(integer sim_meshMe)
!!
!! ARGUMENTS
!!
!!    sim_meshMe      Current Processor Number
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls the RuntimeParameters_get routine for initialization.
!!
!!
!!***

subroutine Simulation_init()
  
  use Simulation_data
  use Driver_interface, ONLY : Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none
#include "Flash.h"
#include "Eos.h"
#include "constants.h"

  

  call Driver_getMype(MESH_COMM, sim_meshMe)

  call RuntimeParameters_get( 'smallx', sim_smallX)
  call RuntimeParameters_get( 'rho_ambient', sim_rhoAmbient)
  call RuntimeParameters_get( 'c_ambient', sim_cAmbient)
  call RuntimeParameters_get( 'mach', sim_mach)
  call RuntimeParameters_get('gamma', sim_gamma)

  sim_vecLen = 1
  sim_mode = MODE_DENS_PRES


end subroutine Simulation_init
