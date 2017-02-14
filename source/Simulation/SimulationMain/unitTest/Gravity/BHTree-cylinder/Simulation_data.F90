!!****if* source/Simulation/SimulationMain/unitTest/Gravity/BHTree-cylinder/Simulation_data
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
  integer, save :: sim_nSubZones, sim_nProfile
  integer, save :: sim_myPE, sim_Comm
  real, save :: sim_Lx, sim_Ly, sim_Lz
  real, save :: sim_boltz, sim_mH, sim_pi, sim_newt
  real, save :: smallp, smlrho, smallX
  real, save :: xmax, xmin, ymax, ymin, zmax, zmin
  real, save :: sim_temp_a, sim_temp_c, sim_press_a, sim_dens_c
  real, save :: abar_1, gamma_1, abar_2, gamma_2

!  real, allocatable, save :: sim_zProf(:), sim_rhoProf(:), sim_pProf(:)
!  real, allocatable, save :: sim_vProf(:), sim_xProf(:), sim_PhiProf(:)

  !! *** Some other module variables for the test  *** !!
  real, save :: sim_absErrMax, sim_relErrMax
  real, save :: sim_solutionErrorTolerance1
  real, save :: sim_solutionErrorTolerance2
  real, save :: sim_corrPot

  !! *** EOS Parameters *** !!
  real, save, dimension(EOS_NUM) :: sim_eosArr
  integer, save :: sim_vecLen, sim_mode
  real, save :: sim_abar_1, sim_gamma_1, sim_abar_2, sim_gamma_2

end module Simulation_data


