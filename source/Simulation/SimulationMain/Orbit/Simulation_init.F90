!!****if* source/Simulation/SimulationMain/Orbit/Simulation_init
!!
!! NAME
!!  Simulation_init
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get routine for initialization.
!!  Initializes initial conditions for Particle Orbit problem.
!!
!! ARGUMENTS
!!
!!   
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY:  PhysicalConstants_get

  implicit none

#include "constants.h"
#include "Flash.h"

  

  call RuntimeParameters_get('xmin',sim_xMin)
  call RuntimeParameters_get('ymin',sim_yMin)
  call RuntimeParameters_get('zmin',sim_zMin)
  call RuntimeParameters_get('xmax',sim_xMax)
  call RuntimeParameters_get('ymax',sim_yMax)
  call RuntimeParameters_get('zmax',sim_zMax)
  call RuntimeParameters_get('separation', sim_Sep)
  call RuntimeParameters_get('ext_field', sim_extField)
  call RuntimeParameters_get('num_particles', sim_nPtot)

  call PhysicalConstants_get('Newton', sim_Newton)

  sim_xPosns = 0.5*(sim_xMin+sim_xMax) + 0.5*sim_Sep*(/-1., 1./)
  sim_yPosns = 0.5*(sim_yMin+sim_yMax)*(/1., 1./)
  sim_zPosns = 0.5*(sim_zMin+sim_zMax)*(/1., 1./)

  sim_EPSX = (sim_xMax-sim_xMin) * EPS
  sim_EPSY = (sim_yMax-sim_yMin) * EPS
  sim_EPSZ = (sim_zMax-sim_zMin) * EPS

  if (sim_extField) call RuntimeParameters_get('ptmass', sim_ptMass)

end subroutine Simulation_init
