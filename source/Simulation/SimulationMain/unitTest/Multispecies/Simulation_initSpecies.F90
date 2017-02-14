!!****if* source/Simulation/SimulationMain/unitTest/Multispecies/Simulation_initSpecies
!!
!! NAME
!!
!!  Simulation_initSpecies
!!
!!
!! SYNOPSIS
!!  Simulation_initSpecies()
!!
!! DESCRIPTION
!!
!!  This routine initializes the species and species values needed for a 
!!  unit test of the Multispecies unit.
!!  This routine is called from Multispecies_init.
!!
!!
!!
!!***

subroutine Simulation_initSpecies()
  use Multispecies_interface, ONLY : Multispecies_setProperty

  implicit none

#include "Flash.h"
#include "Multispecies.h"

!! Set up different arbitrary species, with nice round numbers
!!  to make math easier to check.  Make some of them intentionally
!!  non-physical to check error-handling abilities

!  Similar to SF6
  call Multispecies_setProperty(BIG_SPEC, A, 100.)
  call Multispecies_setProperty(BIG_SPEC, Z, 10.)
  call Multispecies_setProperty(BIG_SPEC, GAMMA, 1.0)
! Similar to Air
  call Multispecies_setProperty(SMAL_SPEC, A, 50.0)
  call Multispecies_setProperty(SMAL_SPEC, Z, 5.0)
  call Multispecies_setProperty(SMAL_SPEC, GAMMA, 0.5)
! Has zero properties to mess with your mind
  call Multispecies_setProperty(ZERO_SPEC, A, 20.0)
  call Multispecies_setProperty(ZERO_SPEC, Z, 2.0)
  call Multispecies_setProperty(ZERO_SPEC, GAMMA, 0.2)
! Has negative properties to mess with your mind
  call Multispecies_setProperty(NEG_SPEC, A, 10.0)
  call Multispecies_setProperty(NEG_SPEC, Z, -10.0)
  call Multispecies_setProperty(NEG_SPEC, GAMMA, -1.0)

end subroutine Simulation_initSpecies
