!!****if* source/Grid/GridParticles/GridParticlesMove/paramesh/VirtualParticles/gr_ptVPData
!!
!! NAME
!!
!!  gr_ptVPData
!!
!! SYNOPSIS
!!
!!  use gr_ptVPData
!!
!! DESCRIPTION
!!
!!  Data module for the variables when using virtual particles
!!  with paramesh
!!
!! ARGUMENTS
!!
!!
!!***

#include "constants.h"
Module gr_ptVPData
  real,dimension(LOW:HIGH,MDIM),save :: gr_ptVPBndBox
  real,dimension(MDIM):: gr_ptVPDeltas
  integer, save :: gr_ptVPMaxCount=100
  real, save :: gr_ptVPBufferFactor=1.
end Module gr_ptVPData
