!!****if* source/Simulation/SimulationMain/IsentropicVortex/Simulation_data
!!
!!  NAME
!!    Simulation_data
!!
!!  SYNOPSIS
!!    use Simulation_data
!!
!!  DESCRIPTION
!!
!!    Stores the local data for Simulation setup: IsentropicVortex
!!
!!***
Module Simulation_data
!
  implicit none

#include "Eos.h"
#include "Flash.h"

  real,    save :: sim_gamma
  real,    save :: sim_uAmbient, sim_vAmbient
  real,    save :: sim_vortexStrength
  real,    save :: sim_xctrTrue, sim_yctrTrue
  integer, save :: sim_nxSubint, sim_nySubint

  real,    save :: sim_imax, sim_jmax, sim_imin, sim_jmin
  real,    save :: sim_imidDomain, sim_jmidDomain, sim_diDomain, sim_djDomain

  real,    save :: sim_tStarAmbient, sim_rbar, sim_constAmbient
  real,    save :: sim_smlrho, sim_smallx
  real, save,dimension(EOS_NUM) :: sim_eosData
  real,save, dimension(NSPECIES) :: sim_eosMassFr
  integer, save :: sim_meshMe
end Module Simulation_data
