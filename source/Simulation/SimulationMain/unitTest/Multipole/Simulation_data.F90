!!****if* source/Simulation/SimulationMain/unitTest/Multipole/Simulation_data
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
!!  Stores the local data for the grid multipole solver unit test, based on
!!  the gravitational MacLaurin spheroid problem.
!!
!!***

module Simulation_data

  implicit none
  
  integer, save :: sim_initGeometry
  integer, save :: sim_globalMe
  integer, save :: sim_nsubzones
  real,    save :: sim_passTolerance              ! maximum percentage error allowed for success
  real,    save :: sim_e                          ! eccentricity
  real,    save :: sim_density
  real,    save :: sim_xctr, sim_yctr, sim_zctr
  real,    save :: sim_Newton, sim_pi
  real,    save :: sim_nsubinv
  real,    save :: sim_a1, sim_a3, sim_a1inv, sim_a3inv
  real,    save :: sim_smallRho

  integer, parameter :: GRID_3DCARTESIAN    =  1, &
                        GRID_3DCYLINDRICAL  =  2, &
                        GRID_2DCYLINDRICAL  =  3, &
                        GRID_1DSPHERICAL    =  4, &
                        GRID_2DSPHERICAL    =  5

end module Simulation_data
