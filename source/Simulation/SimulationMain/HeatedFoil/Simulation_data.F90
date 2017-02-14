!!****if* source/Simulation/SimulationMain/HeatedFoil/Simulation_data
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

  real, save :: sim_foilRadius
  real, save :: sim_foilThickness
  real, save :: sim_foilZPosition
  real, save :: sim_thotFoil
  real, save :: sim_teleRDecayFoil
  real, save :: sim_teleZDecayFoil
  real, save :: sim_rhoFoil
  real, save :: sim_teleFoil 
  real, save :: sim_tionFoil 
  real, save :: sim_tradFoil 
  real, save :: sim_rhoVacu  
  real, save :: sim_teleVacu 
  real, save :: sim_tionVacu 
  real, save :: sim_tradVacu 
  real, save :: sim_smallX


end module Simulation_data


