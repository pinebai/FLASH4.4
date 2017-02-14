!!****if* source/Simulation/SimulationMain/magnetoHD/FieldLoop/Simulation_init
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
!!  none
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!
!!  Reference: Gardiner & Stone JCP 205(2005),509-539
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

  call RuntimeParameters_get('rx',      sim_rx)
  call RuntimeParameters_get('ry',      sim_ry)
  call RuntimeParameters_get('U_initial',   sim_U_initial)
  call RuntimeParameters_get('velz_initial',sim_velz_initial)
  call RuntimeParameters_get('Az_initial',  sim_Az_initial)
  call RuntimeParameters_get('R_fieldLoop', sim_fieldLoopRadius)
  call RuntimeParameters_get('killdivb', sim_killdivb)
  call RuntimeParameters_get('smallp',  sim_smallP)

  sim_gCell = .true.
  
end subroutine Simulation_init
