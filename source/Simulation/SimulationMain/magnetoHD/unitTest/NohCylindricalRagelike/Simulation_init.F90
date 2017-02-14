!!****if* source/Simulation/SimulationMain/magnetoHD/unitTest/NohCylindricalRagelike/Simulation_init
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
!!   
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get routine for initialization.
!!
!!  Reference: Velikovich et al. Phys. Plasmas 19 (2012), 012707
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none

#include "constants.h"
#include "Flash.h"

  

  call RuntimeParameters_get('gamma',   sim_gamma)
  call RuntimeParameters_get('xmin',    sim_xMin)
  call RuntimeParameters_get('ymin',    sim_yMin)
  call RuntimeParameters_get('zmin',    sim_zMin)
  call RuntimeParameters_get('xmax',    sim_xMax)
  call RuntimeParameters_get('ymax',    sim_yMax)
  call RuntimeParameters_get('zmax',    sim_zMax)

  call RuntimeParameters_get('unit_density',     sim_UnitDensity)
  call RuntimeParameters_get('unit_velocity',    sim_UnitVelocity)
  call RuntimeParameters_get('unit_length',      sim_UnitLength)

  call RuntimeParameters_get('killdivb', sim_killdivb)
  call RuntimeParameters_get('smallp',   sim_smallP)

  sim_gCell = .true.
  
end subroutine Simulation_init
