!!****if* source/Simulation/SimulationMain/DustCollapse/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data
!!
!! DESCRIPTION
!!
!!  Store the simulation data for the DustCollapse problem
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!   sim_initRad                  Initial radius of cloud
!!   sim_initDens                 Initial density of cloud
!!   sim_tAmbient                 Initial ambient temperature (everywhere)
!!   sim_iCtr,sim_jCtr, sim_kCtr  Coordinates of the center of the cloud
!!
!!***
Module Simulation_data
  integer,parameter :: N_prof = 1000
  real, save    :: sim_imin, sim_imax, sim_jmin, sim_jmax, sim_kmin, sim_kmax
  real, save    :: sim_smalle, sim_smallp, sim_initRad, sim_gamma, sim_tAmbient
  real, save    :: sim_initDens, sim_smlrho, sim_ictr, sim_jctr, sim_kctr
  real, save    :: sim_presFrac
  real,save, dimension(N_prof):: sim_rProf, sim_rhoProf, sim_pProf,sim_vProf
  integer, save :: sim_meshMe
end Module Simulation_data

