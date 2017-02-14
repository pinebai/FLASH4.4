!!****if* source/Simulation/SimulationMain/magnetoHD/BlastBS/Simulation_init
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
!! ARGUMENTS
!!
!!  none
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use Driver_interface, ONLY : Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none

#include "constants.h"
#include "Flash.h"

  
  call Driver_getMype(MESH_COMM, sim_meshMe)

  call RuntimeParameters_get('gamma',   sim_gamma)
  call RuntimeParameters_get('xmin',    sim_xMin)
  call RuntimeParameters_get('ymin',    sim_yMin)
  call RuntimeParameters_get('zmin',    sim_zMin)
  call RuntimeParameters_get('xmax',    sim_xMax)
  call RuntimeParameters_get('ymax',    sim_yMax)
  call RuntimeParameters_get('zmax',    sim_zMax)
  call RuntimeParameters_get('xCtr',    sim_xCtr)
  call RuntimeParameters_get('yCtr',    sim_yCtr)
  call RuntimeParameters_get('zCtr',    sim_zCtr)
  call RuntimeParameters_get('Radius',  sim_Radius)
#if defined(MAGX_VAR) && defined(MAGY_VAR) && defined(MAGZ_VAR)
  call RuntimeParameters_get('killdivb',sim_killdivb)
  call RuntimeParameters_get('Bx0',     sim_Bx0)
#endif
  call RuntimeParameters_get('smallp',  sim_smallP)
  call RuntimeParameters_get('smallx',  sim_smallX)

  sim_gCell = .true.
  
end subroutine Simulation_init
