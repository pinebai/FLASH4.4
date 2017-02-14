!!****if* source/physics/Eos/EosNuclear/eos_nucData
!!
!! NAME
!!
!!  eos_nucData
!!
!! 
!! SYNOPSIS
!!
!!  use eos_nucData
!!
!! DESCRIPTION
!!
!!  General parameters (non-array) for EOS Nuclear
!!
!! ARGUMENTS
!!
!!
!!*** 

module eos_nucData

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"
#include "Eos_map.h"

  logical, save :: eos_postBounce = .FALSE.
  real, save :: eos_bounceTime = 0.0
  real, save :: eos_centralDens, eos_centralEntr
  integer, save :: eos_nstep

  ! The entropy at which bounce will be flagged:
  real, parameter :: eos_shockEntr = 3.0 
  ! The minimum central density at which bounce will be flagged:
  real, parameter :: eos_bounceDens = 2.0e14

  integer, save :: eos_meshComm
  integer, save :: eos_meshMe

  logical, save :: eos_restart

end module eos_nucData
