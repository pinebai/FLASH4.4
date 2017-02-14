!!****if* source/Simulation/SimulationMain/HeatexchangeIonEle/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!  Use Simulation_data
!!
!! DESCRIPTION
!!
!!  Store the simulation data
!!
!! 
!!***

module Simulation_data

  implicit none

  !! *** Runtime Parameters *** !!

  real, save :: sim_xCenter, sim_yCenter, sim_zCenter
  real, save :: sim_gamma, sim_smallP, sim_smallX
  integer, save ::  sim_orientation
  real, save :: sim_rhoInit, sim_toffset, sim_Q, sim_tempBackground

  logical, save :: sim_gCell

  real, save :: sim_ionTemp, sim_eleTemp, sim_radTemp
  real, save :: sim_CvIon, sim_CvEle
  real, save :: sim_initialCondTemperatureExponent
  real, save :: sim_condTemperatureExponent, sim_alpha, sim_xi0, sim_xfInitial

  real, save :: sim_anaTol, sim_anaSmallT
  ! maximum number of iterations for the Newton loop to find temperate diff from time
  integer, save :: sim_anaMaxNewton
  real,    save :: sim_Avogadro, sim_kBoltzmann, sim_eleCharge, sim_eMassInUAmu
  real, save :: sim_memi, sim_relA, sim_dynamicZ

  integer,save :: sim_schemeOrder
  real,save :: sim_maxTolCoeff0,sim_maxTolCoeff1,sim_maxTolCoeff2,sim_maxTolCoeff3

integer, save :: sim_meshMe
end module Simulation_data


