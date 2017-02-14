!!****if* source/Simulation/SimulationMain/unitTest/Gravity/BHTree-jeans/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!   
!!   use Simulation_data
!!
!! DESCRIPTION
!!
!!  Stores the local data for Simulation setup: BHTree-jeans
!!
!!
!!***


module Simulation_data

  implicit none
#include "Eos.h"

  !! *** Runtime Parameters *** !!
  integer, save :: sim_myPE, sim_Comm
  real, save :: sim_rho0, sim_p0, sim_T0, sim_delta, sim_hx, sim_hy, sim_hz
  real, save :: sim_Lx, sim_Ly, sim_Lz, sim_newton
  real, save :: sim_boltz, sim_mH, sim_pi
  real, save :: smallp, smlrho, smallX
  real, save :: sim_solutionErrorTolerance1
  real, save :: sim_solutionErrorTolerance2

  !! *** Some other module variables for the test  *** !!
  real, save :: sim_absErrMax, sim_relErrMax

  !! *** EOS Parameters *** !!
  real, save, dimension(EOS_NUM) :: sim_eosArr
  integer, save :: sim_vecLen, sim_mode
  real, save :: sim_abar_1, sim_gamma_1

end module Simulation_data


