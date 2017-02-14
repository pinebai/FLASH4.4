!!****if* source/Simulation/SimulationMain/unitTest/Gravity/BHTree/Simulation_data
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
!!  Stores the local data for Simulation setup: StirTurb
!!
!!
!!***


module Simulation_data

  implicit none
#include "Eos.h"

  !! *** Runtime Parameters *** !!
  character(len=255),save :: sim_radprof_file
  integer, save :: sim_pertType, sim_spharm_l1, sim_spharm_m1
  integer, save :: sim_nSubZones, sim_nProfile
  integer, save :: sim_myPE, sim_Comm
  real, save :: sim_xCenter, sim_yCenter, sim_zCenter
  real, save :: sim_vx, sim_vy, sim_vz
  real, save :: sim_Lx, sim_Ly, sim_Lz
  real, save :: sim_boltz, sim_mH, sim_pi
  real, save :: sim_pertamp, sim_velamp
  real, save :: smallp, smlrho, smallX
  real, save :: plgndr_min, plgndr_max, plgndr_norm
  real, save :: jeans_ref, jeans_deref
  real, save :: sim_solutionErrorTolerance1
  real, save :: sim_solutionErrorTolerance2

  !! *** Some other module variables for the test  *** !!
  real, save :: sim_absErrMax, sim_relErrMax

  real, allocatable, save :: sim_rProf(:), sim_rhoProf(:), sim_pProf(:)
  real, allocatable, save :: sim_vProf(:), sim_xProf(:), sim_PhiProf(:)


  !! *** EOS Parameters *** !!
  real, save, dimension(EOS_NUM) :: sim_eosArr
  integer, save :: sim_vecLen, sim_mode
  real, save :: sim_abar_1, sim_gamma_1, sim_abar_2, sim_gamma_2

end module Simulation_data


