!!****if* source/Simulation/SimulationMain/SedovSelfGravity/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data 
!!
!!  DESCRIPTION
!!
!!  Stores the local data for Simulation setup: Sedov
!!  
!! PARAMETERS
!!
!!  p_ambient       Initial ambient pressure
!!  rho_Ambient     Initial ambient density
!!  exp_energy      Explosion energy (distributed over 2^dimen central zones)
!!  t_init          Initial time since explosion
!!  sim_nsubzones      Number of `sub-zones' in cells for applying 1d profile
!!
!!
!!***

Module Simulation_data
#include "Flash.h"
  integer, parameter               :: SIM_NPROFILE = 1000
  integer, save                    :: sim_nsubzones
  real   , save                    :: sim_inSubzones, sim_inSubzm1
  real, save, DIMENSION(NSPECIES)  :: sim_xn

  real, dimension(SIM_NPROFILE), save :: sim_rProf, sim_rhoProf, sim_pProf
  real, dimension(SIM_NPROFILE), save :: sim_vProf

  real, save   :: sim_drProf,sim_gamma,sim_smlrho, sim_smallp
  integer, save :: sim_meshMe
end Module Simulation_data




