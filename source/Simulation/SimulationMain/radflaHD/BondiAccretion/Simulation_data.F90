!!****if* source/Simulation/SimulationMain/radflaHD/BondiAccretion/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data 
!!
!! DESCRIPTION
!!
!!  Stores the simulation data
!!
!!***

module Simulation_data

  implicit none
#include "constants.h"
#include "Flash.h"
  
  integer, save :: nvar_stored
  integer, parameter :: n1d_max = 10000 ! Max number of lines a file can have
  integer, save :: n1d_total ! Actual number of lines, calculated after input
  real, save :: sim_smlrho, sim_smallt,sim_smallx,sim_smallP
  real, save :: sim_smallEele
  character(len=80),save :: model_file
  real,save :: xzn(n1d_max)
  real,save :: model_1d(n1d_max,NUNK_VARS)
  real, save    :: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
  real, save  :: sim_windVel, sim_massLoss, sim_velMult
  integer, save :: nsub
  character (len=4), save :: unklabels(UNK_VARS_BEGIN:UNK_VARS_END)

  integer, save :: sim_meshMe

  logical, save :: sim_restart, sim_burnUpdateEint
  real, save :: sim_pointMass, sim_holeRadius
  real, save :: sim_shelldens, sim_rinner, sim_router, sim_bombRad
  real, save :: sim_ExpEner,sim_coremass,sim_bombRadIn
  real, save :: sim_gamma, sim_TradInitScaleFactor
  real, save :: sim_tele, sim_trad, sim_tion
  real, save :: sim_bondiRadius
  ! Physical constants:
  real, save :: sim_radconst ! Radiation constant
  real, save :: sim_speedlt ! Speed of light
  real, save :: sim_speedlt3 ! a third thereof
  real, save :: sim_t_vac,sim_t_s,sim_rho_vac,sim_rho_s,sim_r_s,sim_steep,sim_rt_s,sim_sb
 
  logical, save :: sim_usePnotT
  logical, save :: sim_plotScaledPressures, sim_useMGD

  real, save :: sim_usedPointMass, sim_luminosity, sim_transOpacConst
  real, save :: sim_soundSpeedInf

  real, save :: sim_holeRad,sim_shelltempfac,sim_accretionRate

  logical, save :: sim_initializeAnalytic

  logical, save :: sim_staticGpot, sim_shellcond,sim_paircond

end module Simulation_data
