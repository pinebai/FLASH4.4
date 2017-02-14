!!****if* source/Simulation/SimulationMain/RHD_Sod/Simulation_init
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
!!  Initializes all the parameters needed for the Sod shock tube
!!  problem
!!
!! ARGUMENTS
!!
!!  
!!
!! PARAMETERS
!!
!!  sim_rhoLeft    Density in the left part of the grid
!!  sim_rhoRight   Density in the right part of the grid
!!  sim_pLeft      Pressure  in the left part of the grid
!!  sim_pRight     Pressure  in the righ part of the grid
!!  sim_uLeft      fluid velocity in the left part of the grid
!!  sim_uRight     fluid velocity in the right part of the grid
!!  sim_xangle     Angle made by diaphragm normal w/x-axis (deg)
!!  sim_ yangle    Angle made by diaphragm normal w/y-axis (deg)
!!  sim_posnR      Point of intersection between the shock plane and the x-axis
!!
!!***

subroutine Simulation_init()
  
  use Simulation_data
  use Driver_interface, ONLY : Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp
  implicit none
#include "Flash.h"
#include "constants.h"

  
  call Driver_getMype(MESH_COMM, sim_meshMe)

  call RuntimeParameters_get('smallp', sim_smallP)
  call RuntimeParameters_get('smallx', sim_smallX) 
  
  call RuntimeParameters_get('gamma', sim_gamma)
  
  call RuntimeParameters_get('sim_rhoLeft',  sim_rhoLeft)
  call RuntimeParameters_get('sim_rhoRight', sim_rhoRight)
  
  call RuntimeParameters_get('sim_pLeft',  sim_pLeft)
  call RuntimeParameters_get('sim_pRight', sim_pRight)
  
  call RuntimeParameters_get('sim_uLeft',  sim_uLeft)
  call RuntimeParameters_get('sim_uRight', sim_uRight)
  call RuntimeParameters_get('sim_vLeft',  sim_vLeft)
  call RuntimeParameters_get('sim_vRight', sim_vRight)
  call RuntimeParameters_get('sim_wLeft',  sim_wLeft)
  call RuntimeParameters_get('sim_wRight', sim_wRight)
  
  call RuntimeParameters_get('sim_xangle', sim_xAngle)
  call RuntimeParameters_get('sim_yangle', sim_yAngle)
  
  call RuntimeParameters_get('sim_posn', sim_posn)


  call Logfile_stamp( "initializing Sod problem",  &
       "[Simulation_init]")
     

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







