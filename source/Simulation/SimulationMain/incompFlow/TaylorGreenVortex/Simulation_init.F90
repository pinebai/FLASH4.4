!!****if* source/Simulation/SimulationMain/incompFlow/TaylorGreenVortex/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_init()
!!
!! ARGUMENTS
!!
!!   none
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!  Initializes initial conditions for INS-isotropic turbulence problem.
!!
!!***

#include "constants.h"

subroutine Simulation_init()

  use Driver_interface, ONLY : Driver_getMype
  use Simulation_data, ONLY : sim_xMin, sim_yMin, &
                              sim_xMax, sim_yMax, &
                              sim_uconv, sim_vconv, sim_meshMe

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

  call Driver_getMype(MESH_COMM, sim_meshMe)

  call RuntimeParameters_get('xmin',    sim_xMin)
  call RuntimeParameters_get('ymin',    sim_yMin)
  call RuntimeParameters_get('xmax',    sim_xMax)
  call RuntimeParameters_get('ymax',    sim_yMax)

  call RuntimeParameters_get('uconv',   sim_uconv)
  call RuntimeParameters_get('vconv',   sim_vconv)

end subroutine Simulation_init
