!!****if* source/Simulation/SimulationMain/radflaHD/CriticalRadShock/Simulation_init
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
  use Logfile_interface, ONLY : Logfile_stamp
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
#include "Flash.h"

  call RuntimeParameters_get('sim_rho' , sim_rho)
  call RuntimeParameters_get('sim_tgas', sim_tgas)
  call RuntimeParameters_get('sim_trad', sim_trad)
  call RuntimeParameters_get('sim_velx', sim_velx)
  call RuntimeParameters_get('xmin', sim_xmin)
  call RuntimeParameters_get('xmax', sim_xmax)
  call RuntimeParameters_get('ymin', sim_ymin)
  call RuntimeParameters_get('ymax', sim_ymax)
  call RuntimeParameters_get('zmin', sim_zmin)
  call RuntimeParameters_get('zmax', sim_zmax)

  call RuntimeParameters_get("gamma",sim_gamma)

end subroutine Simulation_init
