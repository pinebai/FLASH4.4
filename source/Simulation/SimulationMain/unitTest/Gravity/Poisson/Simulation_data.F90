!!****if* source/Simulation/SimulationMain/unitTest/Gravity/Poisson/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!  module Simulation_data()
!!
!! DESCRIPTION
!!
!!  Store the simulation data for unitTesting of Poisson solver
!!   
!! ARGUMENTS
!!
!! PARAMETERS   
!!      Described in the Config file
!!
!! NOTES
!!
!!  No arguments.  All data passed by "use Simulation_data"
!!
!!***

module Simulation_data

  implicit none
#include "constants.h"
#include "Flash.h"
  
  integer, save :: sim_meshMe
  real, save :: sim_smlrho
  integer, save :: sim_subSample

 
  integer, parameter :: Npeak = 12
  real, save :: sim_xctr(Npeak) = & 
       &          (/  -1., -1., 0.28125, 0.5, 0.3046875, 0.3046875, & 
       &             0.375, 0.5625, -0.5, -0.125, 0.296875, 0.5234375 /)
  real, save :: sim_yctr(Npeak) = & 
       &          (/  0.09375, 1., 0.53125, 0.53125, 0.1875, 0.125, & 
       &             0.15625, -0.125, -0.703125, -0.703125, -0.609375, & 
       &             -0.78125 /)
  real, save :: sim_zctr(Npeak) = & 
       &          (/ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0./)
  real, save :: sim_sigma(Npeak) = & 
       &          (/ 4000., 20000., 80000., 16., 360000., 400000., & 
       &             2000., 18200., 128., 49000., 37000., 18900. /)
  

end module Simulation_data
