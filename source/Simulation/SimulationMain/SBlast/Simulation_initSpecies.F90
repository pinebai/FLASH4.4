!!****if* source/Simulation/SimulationMain/SBlast/Simulation_initSpecies
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

  use Simulation_data, ONLY: &
       sim_AIn, sim_ZIn, sim_gammaIn, &
       sim_A1 , sim_Z1 , sim_gamma1 , &
       sim_A2 , sim_Z2 , sim_gamma2
  use Multispecies_interface, ONLY : Multispecies_setProperty
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

#include "Multispecies.h"
#include "Flash.h"


  call RuntimeParameters_get('sim_AIn', sim_AIn)
  call RuntimeParameters_get('sim_ZIn', sim_ZIn)
  call RuntimeParameters_get('sim_gammaIn', sim_gammaIn)

  call RuntimeParameters_get('sim_A1', sim_A1)
  call RuntimeParameters_get('sim_Z1', sim_Z1)
  call RuntimeParameters_get('sim_gamma1', sim_gamma1)

  call RuntimeParameters_get('sim_A2', sim_A2)
  call RuntimeParameters_get('sim_Z2', sim_Z2)
  call RuntimeParameters_get('sim_gamma2', sim_gamma2)

  call Multispecies_setProperty(FLD1_SPEC, A, sim_AIn)
  call Multispecies_setProperty(FLD1_SPEC, Z, sim_ZIn)
  call Multispecies_setProperty(FLD1_SPEC, GAMMA, sim_gammaIn)

  call Multispecies_setProperty(FLD2_SPEC, A, sim_A1)
  call Multispecies_setProperty(FLD2_SPEC, Z, sim_Z1)
  call Multispecies_setProperty(FLD2_SPEC, GAMMA, sim_gamma1)

  call Multispecies_setProperty(FLD3_SPEC, A, sim_A2)
  call Multispecies_setProperty(FLD3_SPEC, Z, sim_Z2)
  call Multispecies_setProperty(FLD3_SPEC, GAMMA, sim_gamma2)


end subroutine Simulation_initSpecies

