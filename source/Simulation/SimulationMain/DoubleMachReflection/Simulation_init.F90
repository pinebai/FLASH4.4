!!****if* source/Simulation/SimulationMain/DoubleMachReflection/Simulation_init
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
!!  sim_uLeft      fluid x-velocity in the left part of the grid
!!  sim_uRight     fluid x-velocity in the right part of the grid
!!  sim_vLeft      fluid y-velocity in the left part of the grid
!!  sim_vRight     fluid y-velocity in the right part of the grid
!!  sim_xangle     Angle made by diaphragm normal w/x-axis (deg)
!!  sim_posn       Point of intersection between the shock plane and the x-axis
!!
!!***

subroutine Simulation_init()
  
  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp
  use Driver_interface, ONLY : Driver_getMype
  implicit none
#include "Flash.h"
#include "constants.h"
  

  call RuntimeParameters_get('smallp', sim_smallP)
  call RuntimeParameters_get('smallx', sim_smallX) 
  
  call RuntimeParameters_get('gamma', sim_gamma)
  
  call RuntimeParameters_get('sim_rhoLeft', sim_rhoLeft)
  call RuntimeParameters_get('sim_rhoRight', sim_rhoRight)
  
  call RuntimeParameters_get('sim_pLeft', sim_pLeft)
  call RuntimeParameters_get('sim_pRight', sim_pRight)
  
  call RuntimeParameters_get('sim_uLeft', sim_uLeft)
  call RuntimeParameters_get('sim_uRight', sim_uRight)

  call RuntimeParameters_get('sim_vLeft', sim_vLeft)
  call RuntimeParameters_get('sim_vRight', sim_vRight)
  
  call RuntimeParameters_get('sim_xangle', sim_xAngle)
  
  call RuntimeParameters_get('sim_posn', sim_posn)

  call RuntimeParameters_get('xmin',      sim_xmin)
  call RuntimeParameters_get('xmax',      sim_xmax)
  call RuntimeParameters_get('ymin',      sim_ymin)
  call RuntimeParameters_get('ymax',      sim_ymax)
  call Driver_getMype(MESH_COMM, sim_meshMe)


  call Logfile_stamp( "initializing Sod problem",  &
       "[Simulation_init]")
     

  ! convert the shock angle paramters
  sim_xAngle = sim_xAngle * 0.0174532925 ! Convert to radians.


end subroutine Simulation_init







