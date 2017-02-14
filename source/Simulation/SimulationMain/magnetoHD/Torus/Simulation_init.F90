!!****if* source/Simulation/SimulationMain/magnetoHD/Torus/Simulation_init
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

  call RuntimeParameters_get('R_min',     sim_rMin)
  call RuntimeParameters_get('R_max',     sim_rMax)
  call RuntimeParameters_get('den_max',   sim_denMax)
  call RuntimeParameters_get('den_cut',   sim_denCut)
  call RuntimeParameters_get('BETA',      sim_beta)
  call RuntimeParameters_get('D_Con',     sim_dCon)
  call RuntimeParameters_get('T_Con',     sim_tCon)
  call RuntimeParameters_get('R_0',       sim_r0)
  call RuntimeParameters_get('R_Sphere',  sim_rSphere)
  
  call RuntimeParameters_get('gamma',   sim_gamma)
  call RuntimeParameters_get('xmin',    sim_xMin)
  call RuntimeParameters_get('ymin',    sim_yMin)
  call RuntimeParameters_get('zmin',    sim_zMin)
  call RuntimeParameters_get('xmax',    sim_xMax)
  call RuntimeParameters_get('ymax',    sim_yMax)
  call RuntimeParameters_get('zmax',    sim_zMax)


  call RuntimeParameters_get('killdivb', sim_killdivb)
  call RuntimeParameters_get('smallp',  sim_smallP)

  sim_gCell = .true.
  
end subroutine Simulation_init
