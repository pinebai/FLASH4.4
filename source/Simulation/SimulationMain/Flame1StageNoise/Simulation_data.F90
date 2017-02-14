module Simulation_data
#include "Flash.h"
#include "Eos.h"
  real, save :: sim_rhoAmbient, sim_tAmbient
  logical, save :: sim_ignite, sim_pseudo1d
  real, save :: sim_fracPerturb
  real, save :: sim_xctrPerturb, sim_yctrPerturb, sim_zctrPerturb, sim_theta

  real, save :: sim_xmin, sim_xmax, sim_ymin, sim_ymax

  real,    save           :: sim_laminarWidth, sim_flamespeed
  real, dimension(EOS_NUM):: sim_eosData_u, sim_eosData_b

end module Simulation_data
