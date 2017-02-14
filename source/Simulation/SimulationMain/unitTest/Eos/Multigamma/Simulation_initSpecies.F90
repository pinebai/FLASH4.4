!!****if* source/Simulation/SimulationMain/unitTest/Eos/Multigamma/Simulation_initSpecies
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
!!  This routine will initialize the simulation with the species 
!!  in the Config files. This particular implementation has air
!!  and SF6, and they have two different gamma values for the 
!!  ideal gas gamma law. 
!!
!!
!!
!!***

subroutine Simulation_initSpecies()
  use Multispecies_interface, ONLY : Multispecies_setProperty

implicit none
#include "Flash.h"
#include "Multispecies.h"


  call Multispecies_setProperty(HE4_SPEC, A, 4.)
  call Multispecies_setProperty(HE4_SPEC, Z, 2.)
  call Multispecies_setProperty(HE4_SPEC, GAMMA, 1.66)

  call Multispecies_setProperty(NI56_SPEC, A, 56.)
  call Multispecies_setProperty(NI56_SPEC, Z, 28.)
  call Multispecies_setProperty(NI56_SPEC, GAMMA, 1.4)

end subroutine Simulation_initSpecies
