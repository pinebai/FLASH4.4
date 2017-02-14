!!****if* source/Simulation/SimulationMain/unitTest/Gravity/BHTree-layer/Simulation_data
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
!!  Stores the local data for Simulation setup: BHTree-layer
!!
!!
!!***


module Simulation_data

  implicit none
#include "Eos.h"

  !! *** Runtime Parameters *** !!
  character(len=255),save :: sim_prof_file
  integer, save :: sim_nSubZones, sim_nProfile, sim_dir
  integer, save :: sim_myPE, sim_Comm
  real, save :: sim_zMidplane
  real, save :: sim_Lx, sim_Ly, sim_Lz
  real, save :: sim_boltz, sim_mH, sim_pi
  real, save :: smallp, smlrho, smallX
  real, save :: sim_solutionErrorTolerance1
  real, save :: sim_solutionErrorTolerance2

  real, allocatable, save :: sim_zProf(:), sim_rhoProf(:), sim_pProf(:)
  real, allocatable, save :: sim_vProf(:), sim_xProf(:), sim_PhiProf(:)

  !! *** Some other module variables for the test  *** !!
  real, save :: sim_absErrMax, sim_relErrMax

  !! *** EOS Parameters *** !!
  real, save, dimension(EOS_NUM) :: sim_eosArr
  integer, save :: sim_vecLen, sim_mode
  real, save :: sim_abar_1, sim_gamma_1

end module Simulation_data


