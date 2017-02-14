!!****if* source/Simulation/SimulationMain/Blast2/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!  
!! use Simulation_data
!! 
!! DESCRIPTION
!!  Stores Simulation unit scope data
!!  for applying initial conditions to
!!  to the Blast2 problem. 
!!   
!!
!!***

module Simulation_data

  implicit none

  !! *** Runtime Parameters *** !!
  real, save :: sim_gamma, sim_smallP, sim_smallX, sim_smallE, sim_smallRho
  real, save :: sim_rhoLeft, sim_rhoMid, sim_rhoRight, & 
                sim_pLeft, sim_pMid, sim_pRight, & 
                sim_uLeft, sim_uMid, sim_uRight, & 
                sim_xAngle, sim_yAngle, sim_posnL, sim_posnR

  !! *** Variables pertaining to this Simulation *** !!
  real, save :: sim_xCos, sim_yCos, sim_zCos

  real, save :: sim_eMassInUAmu
integer, save :: sim_meshMe
end module Simulation_data


