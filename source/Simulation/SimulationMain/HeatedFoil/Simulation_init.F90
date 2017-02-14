!!****if* source/Simulation/SimulationMain/HeatedFoil/Simulation_init
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
!!  Initializes all the parameters needed for a particular simulation
!!
!!
!! ARGUMENTS
!!
!!  
!!
!! PARAMETERS
!!
!!***

subroutine Simulation_init()
  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  
  implicit none

#include "constants.h"

  character(len=MAX_STRING_LENGTH) :: str

  call RuntimeParameters_get('sim_foilRadius', sim_foilRadius)
  call RuntimeParameters_get('sim_foilThickness', sim_foilThickness)
  call RuntimeParameters_get('sim_foilZPosition', sim_foilZPosition)
  call RuntimeParameters_get('sim_thotFoil', sim_thotFoil)
  call RuntimeParameters_get('sim_teleRDecayFoil', sim_teleRDecayFoil)
  call RuntimeParameters_get('sim_teleZDecayFoil', sim_teleZDecayFoil)  
  call RuntimeParameters_get('sim_rhoFoil', sim_rhoFoil)
  call RuntimeParameters_get('sim_teleFoil', sim_teleFoil)
  call RuntimeParameters_get('sim_tionFoil', sim_tionFoil)
  call RuntimeParameters_get('sim_tradFoil', sim_tradFoil)
  call RuntimeParameters_get('sim_rhoVacu', sim_rhoVacu)
  call RuntimeParameters_get('sim_teleVacu', sim_teleVacu)
  call RuntimeParameters_get('sim_tionVacu', sim_tionVacu)
  call RuntimeParameters_get('sim_tradVacu', sim_tradVacu)
  call RuntimeParameters_get('smallX', sim_smallX)

end subroutine Simulation_init
