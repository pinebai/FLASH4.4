!!****if* source/Simulation/SimulationMain/NeiTest/Simulation_init
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
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!  Initializes initial conditions for Sedov Spherical Explosion 
!!  problem.
!!
!! ARGUMENTS
!!
!!   
!!
!! PARAMETERS
!!  sim_rhoAmbient     Initial ambient density
!!  sim_velInit        Initial velocity
!!  sim_tAmbient       Initial ambient temperature
!!  sim_tPerturb       Perturbation temperature
!!  sim_radius         Radial position of inner edge of grid (for 1D )
!!
!!
!!***

!There is duplication between this and source/physics/sourceTerms/Ionize/IonizeMain/Ionize_init

subroutine Simulation_init()

  use Simulation_data 
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

  
  
  call RuntimeParameters_get('smallx', sim_smallx)
  call RuntimeParameters_get('vel_init', sim_velInit)
  call RuntimeParameters_get('rho_ambient', sim_rhoAmbient)
  call RuntimeParameters_get('t_ambient', sim_tAmbient)
  call RuntimeParameters_get('t_perturb', sim_tPerturb)
  call RuntimeParameters_get('radius', sim_radius)
  call RuntimeParameters_get('xstep', sim_xstep)


  !Grab some Config file parameters
  call RuntimeParameters_get("vel_init", vel_init)
  
  return
end subroutine Simulation_init
