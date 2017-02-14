!!****if* source/Simulation/SimulationMain/unitTest/ParticlesAdvance/HomologousPassive/Simulation_init
!!
!! NAME
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!  Simulation_init( )
!!
!! DESCRIPTION   
!!   Initialize all the runtime parameters needed for the Particle unitTest
!!
!! ARGUMENTS
!!
!!   
!!
!! PARAMETERS
!!
!!   sim_rho_amb    Gas Density:  Entire domain receives this ambient parameter
!!   sim_p_amb      Gas Pressure: Entire domain receives this ambient parameter
!!   sim_vx_amb     Gas x-velocity:  Dominant flow velocity throughout domain 
!!   sim_vx_multiplier   Half of the domain in y has x-velocity multiplied by this value
!!   sim_seed   Random number seed -- NOT USED please ignore
!!   sim_vx_pert   Scales [-1,1] random number in x direction: set to zero for uniform flow
!!   sim_vy_pert   Scales [-1,1] random number in y direction: set to zero for uniform flow
!!   sim_vz_pert   Scales [-1,1] random number in z direction: set to zero for uniform flow
!!
!!   sim_analyticParticlePositions 
!!   sim_maxTolCoeff0,sim_maxTolCoeff1,sim_maxTolCoeff2,sim_maxTolCoeff3
!!***

subroutine Simulation_init()

  use Simulation_data
  use pt_utData
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none

#include "constants.h"
#include "Flash.h"
  
  
  integer :: i, j, status
 

!! ambient values
  call RuntimeParameters_get("sim_p_amb",   sim_p_amb )
  call RuntimeParameters_get("sim_rho_amb", sim_rho_amb )
  call RuntimeParameters_get("sim_vx_amb",  sim_vx_amb )
  call RuntimeParameters_get("sim_vx_multiplier",  sim_vx_multiplier )

!! velocity perturbations
  call RuntimeParameters_get("sim_seed",   sim_seed )
  call RuntimeParameters_get("sim_vx_pert",   sim_vx_pert )
  call RuntimeParameters_get("sim_vy_pert",   sim_vy_pert )
  call RuntimeParameters_get("sim_vz_pert",   sim_vz_pert )

!! Components of homologous expansion factor: a(t) = a0 + a1*t
  call RuntimeParameters_get("sim_a0",   sim_a0 )
  call RuntimeParameters_get("sim_a1",   sim_a1 )
  call RuntimeParameters_get("sim_analyticParticlePositions",   sim_analyticParticlePositions )
  call RuntimeParameters_get("sim_fakeMapMeshToParticles",      sim_fakeMapMeshToParticles)

  pt_utA0 = sim_a0
  pt_utA1 = sim_a1
  pt_utAnalyticParticlePositions = sim_analyticParticlePositions
  
  call RuntimeParameters_get("tinitial",pt_utInitialSimTime)
  call RuntimeParameters_get("sim_schemeOrder", sim_schemeOrder )
  call RuntimeParameters_get("sim_maxTolCoeff0", sim_maxTolCoeff0 )
  call RuntimeParameters_get("sim_maxTolCoeff1", sim_maxTolCoeff1 )
  call RuntimeParameters_get("sim_maxTolCoeff2", sim_maxTolCoeff2 )
  call RuntimeParameters_get("sim_maxTolCoeff3", sim_maxTolCoeff3 )
  return
end subroutine Simulation_init
