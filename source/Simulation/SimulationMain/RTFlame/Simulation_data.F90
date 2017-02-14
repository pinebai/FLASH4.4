! Dean Townsley 2009

module Simulation_data

  implicit none
#include "Eos.h"

  ! place to store initial data
  real, save, allocatable, dimension(:)  :: sim_dens_i, sim_temp_i, sim_eint_i, sim_ye_i, sim_sumy_i, sim_flam_i
  real, save  :: sim_grav

  real, save  :: sim_last_burned_mass

  integer,save :: HSE_FORWARD = 1
  integer,save :: HSE_BACKWARD = 2
  integer,save :: HSE_CONSTENTR = 3
  integer,save :: HSE_CONSTTEMP = 4
  integer,save :: HSE_SETTEMP = 5

  real, save   :: sim_temp_u, sim_dens_u
  real, save   :: sim_ye_u, sim_sumy_u
  real, save   :: sim_temp_b, sim_dens_b
  real, save   :: sim_ye_b, sim_sumy_b

  real, save   :: sim_eosData_u(EOS_NUM), sim_eosData_b(EOS_NUM)

  real, save   :: sim_x0

  real, save   :: sim_vel_pert_amp, sim_vel_pert_wavelength1

  real, save   :: sim_spert_ampl1, sim_spert_wl1, sim_spert_phase1
  real, save   :: sim_spert_ampl2, sim_spert_wl2, sim_spert_phase2

  logical, save:: sim_refine_uniform_region
  real, save   :: sim_refine_region_size, sim_refine_lead, sim_refine_buf
  real, save   :: sim_refine_region_stepdown_size

  logical, save:: sim_ParticleRefineRegion
  integer, save:: sim_ParticleRefineRegionLevel
  real, save   :: sim_ParticleRefineRegionBottom, sim_ParticleRefineRegionTop

end module Simulation_data
