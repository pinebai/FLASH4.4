!!****if* source/Simulation/SimulationMain/SBlast/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data
!!
!!***

module Simulation_data

  implicit none

  !! *** Runtime Parameters *** !!

  integer :: sim_geo, sim_atmos1, sim_atmos2

  logical :: sim_useE, sim_expIn, sim_ibound

  real :: sim_rhoIn, sim_pIn, sim_EIn, sim_AIn, sim_ZIn, sim_gammaIn, sim_xcIn, sim_ycIn, sim_zcIn, sim_rIn
  real :: sim_rho1, sim_p1, sim_A1, sim_Z1, sim_gamma1, sim_h1, sim_sh1
  real :: sim_rho2, sim_p2, sim_A2, sim_Z2, sim_gamma2, sim_sh2

  real :: xmin, xmax, ymin, ymax, zmin, zmax

integer, save :: sim_meshMe
end module Simulation_data


