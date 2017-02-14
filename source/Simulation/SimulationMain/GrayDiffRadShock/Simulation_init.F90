!!****if* source/Simulation/SimulationMain/GrayDiffRadShock/Simulation_init
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
  call RuntimeParameters_get('sim_temp', sim_temp)
  call RuntimeParameters_get('sim_M0'  , sim_M0)
  call RuntimeParameters_get('sim_P0'  , sim_P0)


end subroutine Simulation_init
