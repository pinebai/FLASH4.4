!!****if* source/Simulation/SimulationMain/MGDStep/Simulation_initSpecies
!!
!! NAME
!!
!!  Simulation_initSpecies
!!
!! SYNOPSIS
!!
!!  Simulation_initSpecies()
!!
!! DESCRIPTION
!!
!!  This routine will initialize the species and species values needed for a 
!!  given setup.   The user should add the 
!!  implementation of this routine to the setups directory of a simulation 
!!  that needs to use the multispecies capabilities of the code.
!!
!!  There two general purpose implementations available in the code, one which sets standard  
!!  isotope properties for the nuclear burning source terms, and another one for the 
!!  Ionization source term.
!!
!!  This routine is called from Multispecies_init, and is called BEFORE
!!  the call to Simulation_init.  
!!
!! SEE ALSO
!!  Multispecies_init
!!  Simulation/SimulationComposition/Simulation_initSpecies
!!
!!***

subroutine Simulation_initSpecies()
  use Multispecies_interface, ONLY : Multispecies_setProperty
  implicit none

#include "Flash.h"
#include "Multispecies.h"

  real :: aelems(MS_MAXELEMS)
  integer :: zelems(MS_MAXELEMS)
  real :: fractions(MS_MAXELEMS)

  call Multispecies_setProperty(HE_SPEC, A, 4.0026032497)
  call Multispecies_setProperty(HE_SPEC, Z, 2.0)
  call Multispecies_setProperty(HE_SPEC, GAMMA, 1.4)
  call Multispecies_setProperty(HE_SPEC, MS_NUMELEMS, 1)

  aelems(1) = 4.0026032497
  call Multispecies_setProperty(HE_SPEC, MS_AELEMS, aelems)

  zelems(1) = 2
  call Multispecies_setProperty(HE_SPEC, MS_ZELEMS, zelems)

  fractions(1) = 1.0
  call Multispecies_setProperty(HE_SPEC, MS_FRACTIONS, fractions)

end subroutine Simulation_initSpecies
