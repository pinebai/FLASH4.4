!!****if* source/Simulation/SimulationMain/FlatPlate/Simulation_init
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
!!  Initializes all the parameters needed for the wind tunnel with a step problem
!!
!! ARGUMENTS
!!
!!  
!!
!! PARAMETERS
!!
!!  sim_pAmbient    Initial ambient pressure
!!  sim_rhoAmbient  Initial ambient density
!!  sim_windVel     Inflow velocity (parallel to x-axis)
!!  gamma           the Gamma EOS thing
!!  smallp          minimum for pressure
!!  smallx          minimum for abundance
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp

  implicit none

#include "constants.h"
#include "Flash.h"

  

  call RuntimeParameters_get('sim_pAmbient', sim_pAmbient)
  call RuntimeParameters_get('sim_rhoAmbient', sim_rhoAmbient)
  call RuntimeParameters_get('sim_rhoBulk', sim_rhoBulk)
  call RuntimeParameters_get('sim_windVelx', sim_windVelx)
  call RuntimeParameters_get('sim_windVely', sim_windVely)
  call RuntimeParameters_get('sim_windVelz', sim_windVelz)
  call RuntimeParameters_get('gamma', sim_gamma)
  call RuntimeParameters_get('smallx', sim_smallX)
  call RuntimeParameters_get('smallp', sim_smallP)
  call RuntimeParameters_get('sim_radius', sim_radius)
  call RuntimeParameters_get('sim_xCtr', sim_xCtr)
  call RuntimeParameters_get('sim_yCtr', sim_yCtr)
  call RuntimeParameters_get('sim_zCtr', sim_zCtr)
  call RuntimeParameters_get('sim_number', sim_number)
  call RuntimeParameters_get('sim_Mach',   sim_Mach)
  call RuntimeParameters_get('sim_xAngle', sim_xAngle)
  call RuntimeParameters_get('sim_yAngle', sim_yAngle)
  call RuntimeParameters_get('sim_zAngle', sim_zAngle)

  sim_gCell = .true.

  call Logfile_stamp("initializing for windtunnel + step", 'run_init')
  write (*,*) "flash:  initializing for wind tunnel + step"

  ! convert the shock angle paramters
  sim_xAngle = sim_xAngle * 0.0174532925 ! Convert to radians.
  sim_yAngle = sim_yAngle * 0.0174532925

  sim_xCos = cos(sim_xAngle)
  
  if (NDIM == 1) then
     sim_xCos = 1.
     sim_yCos = 0.
     sim_zCos = 0.
     
  elseif (NDIM == 2) then
     sim_yCos = sqrt(1. - sim_xCos*sim_xCos)
     sim_zCos = 0.
     
  elseif (NDIM == 3) then
     sim_yCos = cos(sim_yAngle)
     sim_zCos = sqrt( max(0., 1. - sim_xCos*sim_xCos - sim_yCos*sim_yCos) )
  endif

end subroutine Simulation_init
