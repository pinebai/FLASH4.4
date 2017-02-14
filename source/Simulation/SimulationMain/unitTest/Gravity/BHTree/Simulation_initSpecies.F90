!!****if* source/Simulation/SimulationMain/unitTest/Gravity/BHTree/Simulation_initSpecies
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
!!  This routine will initialize the species and species values needed for
!!  the TwoGamma setup, which advects two fluids with different Gamma values
!!
!!***

subroutine Simulation_initSpecies()

  use Simulation_data
  use Multispecies_interface, ONLY : Multispecies_setProperty
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none
!#include "RuntimeParameters.h"
#include "Multispecies.h"
#include "Flash.h"
!#include "Multispecies_interface.h"
  integer :: i


  call RuntimeParameters_get('abar_1', sim_abar_1)
  call RuntimeParameters_get('gamma_1', sim_gamma_1)

  call RuntimeParameters_get('abar_2', sim_abar_2)
  call RuntimeParameters_get('gamma_2', sim_gamma_2)


  call Multispecies_setProperty(FLD1_SPEC, A, sim_abar_1)
  call Multispecies_setProperty(FLD1_SPEC, Z, sim_abar_1)
  call Multispecies_setProperty(FLD1_SPEC, GAMMA, sim_gamma_1)

  call Multispecies_setProperty(FLD2_SPEC, A, sim_abar_2)
  call Multispecies_setProperty(FLD2_SPEC, Z, sim_abar_2)
  call Multispecies_setProperty(FLD2_SPEC, GAMMA, sim_gamma_2)



end subroutine Simulation_initSpecies

