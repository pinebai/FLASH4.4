!
! Dean Townsley 2008
!

module Simulation_data
#include "Flash.h"
#include "constants.h"
#include "Eos.h"
  real, save :: sim_rhoAmbient, sim_tAmbient
  logical, save :: sim_ignite, sim_restartVortex
  logical, save :: sim_inflowVortex, sim_variableInflow
  real, save :: sim_fracPerturb
  real, save :: sim_cFrac, sim_neFrac
  real, save :: sim_sigT, sim_sigP, sim_sigVx, sim_sigVy, sim_sigVz
  real, save :: sim_inflowVx, sim_crossArea
  real, save :: sim_yctrVortex, sim_vortexStrength, sim_vortexSize
  real, save :: sim_xBeginVortex, sim_xEndVortex

  real, save :: sim_vrms
  integer, save :: sim_smooth_level
  character(len=4096), save :: sim_turbfield_filename
  real, dimension(IAXIS:KAXIS,LOW:HIGH), save :: sim_turbfield_bbox

  real, save :: sim_xmin, sim_xmax, sim_ymin, sim_ymax, sim_zmin, sim_zmax

  real, save :: sim_laminarWidth, sim_flamespeed
  real, dimension(EOS_NUM), save :: sim_eosData_u, sim_eosData_b

  real, save :: sim_last_burned_mass


end module Simulation_data
