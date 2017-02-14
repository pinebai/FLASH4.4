!!****if* source/Simulation/SimulationMain/PreHeatPlug/Simulation_init
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
!!  2014/11/10  Add the material for filled gas
!!  2015/4/29   Add the material for LEH window and Washer (named as Wash)
!!              WashThickness > TargetThickness and WashRadius > TargetRadius are needed
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
  
  implicit none

#include "constants.h"
#include "Flash.h"

  real :: xmin, xmax, ymin, ymax
  integer :: lrefine_max, nblockx, nblocky
  character(len=MAX_STRING_LENGTH) :: str

  call RuntimeParameters_get('sim_targetRadius', sim_targetRadius)
  call RuntimeParameters_get('sim_targetHeight', sim_targetHeight)
  call RuntimeParameters_get('sim_vacuumHeight', sim_vacuumHeight)
  call RuntimeParameters_get('sim_targetThickness', sim_targetThickness)        !2014/11/10, by Po-Yu
  call RuntimeParameters_get('sim_gasRadius', sim_gasRadius)                    !2014/11/10, by Po-Yu
  call RuntimeParameters_get('sim_windowsThickness', sim_windowsThickness)      !2014/11/10, by Po-Yu
  call RuntimeParameters_get('sim_windowsRadius', sim_windowsRadius)            !2014/11/10, by Po-Yu  
  call RuntimeParameters_get('sim_washThickness', sim_washThickness)      !2015/4/30, by Po-Yu
  call RuntimeParameters_get('sim_washRadius', sim_washRadius)            !2015/4/30, by Po-Yu    
  call RuntimeParameters_get('sim_plugThickness', sim_plugThickness)            !2015/9/12, by Po-Yu    
  
  call RuntimeParameters_get('sim_rhoTarg', sim_rhoTarg)
  call RuntimeParameters_get('sim_teleTarg', sim_teleTarg)
  call RuntimeParameters_get('sim_tionTarg', sim_tionTarg)
  call RuntimeParameters_get('sim_tradTarg', sim_tradTarg)

  call RuntimeParameters_get('sim_rhoCham', sim_rhoCham)
  call RuntimeParameters_get('sim_teleCham', sim_teleCham)
  call RuntimeParameters_get('sim_tionCham', sim_tionCham)
  call RuntimeParameters_get('sim_tradCham', sim_tradCham)
  
  call RuntimeParameters_get('sim_rhoGas', sim_rhoGas)                  !2014/11/10, by Po-Yu
  call RuntimeParameters_get('sim_teleGas', sim_teleGas)                !2014/11/10, by Po-Yu
  call RuntimeParameters_get('sim_tionGas', sim_tionGas)                !2014/11/10, by Po-Yu
  call RuntimeParameters_get('sim_tradGas', sim_tradGas)                !2014/11/10, by Po-Yu

  call RuntimeParameters_get('sim_rhoLEH', sim_rhoLEH)                  !2015/4/29, by Po-Yu
  call RuntimeParameters_get('sim_teleLEH', sim_teleLEH)                !2015/4/29, by Po-Yu
  call RuntimeParameters_get('sim_tionLEH', sim_tionLEH)                !2015/4/29, by Po-Yu
  call RuntimeParameters_get('sim_tradLEH', sim_tradLEH)                !2015/4/29, by Po-Yu
  
  call RuntimeParameters_get('sim_rhoWash', sim_rhoWash)                  !2015/4/29, by Po-Yu
  call RuntimeParameters_get('sim_teleWash', sim_teleWash)                !2015/4/29, by Po-Yu
  call RuntimeParameters_get('sim_tionWash', sim_tionWash)                !2015/4/29, by Po-Yu
  call RuntimeParameters_get('sim_tradWash', sim_tradWash)                !2015/4/29, by Po-Yu  

  call RuntimeParameters_get('smallX', sim_smallX)

  call RuntimeParameters_get('sim_initGeom', sim_initGeom)

end subroutine Simulation_init
