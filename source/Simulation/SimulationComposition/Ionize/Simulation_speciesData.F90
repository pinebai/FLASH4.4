!!****if* source/Simulation/SimulationComposition/Ionize/Simulation_speciesData
!!
!! NAME
!!
!!  Simulation_speciesData
!!
!! SYNOPSIS
!!
!!  use Simulation_speciesData 
!!
!!  DESCRIPTION
!!
!!  Stores the local data for Simulation setup: Neitest
!!  
!! PARAMETERS
!!
!!  sim_rhoAmbient     Initial ambient density
!!  sim_velInit        Initial velocity
!!  sim_tAmbient       Initial ambient temperature
!!  sim_tPerturb       Perturbation temperature
!!  sim_radius         Radial position of inner edge of grid (for 1D )
!!
!!
!!***

module Simulation_speciesData

#include "Flash.h"
#include "Ionize.h"

  implicit none
  
  !It is assumed that all simulations involving the Ionize unit 
  !have an electron and Hydrogen species.
  integer, parameter :: sim_ELEC = ELEC_SPEC - SPECIES_BEGIN + 1
  integer, parameter :: sim_HYD = H_SPEC - SPECIES_BEGIN + 1


  ! save the parameters that describe this initialization
  integer, parameter :: SPEC_NUM=ION_NELEM, SPEC_NIMAX=ION_NIMAX

  !DEV : The appropriate type is commented out.  This is because our only 
  !test for correctness is to validate the code against FLASH2 implementation 
  !and to generate the same Neitest graphs as in FLASH2 userguide. 
  !For consistency, the effective value of sim_specEpRatio is now 1.0,  
  !due to it being cast to integer (as in FLASH2).
  !  real, parameter :: sim_specEpRatio=1.213
  integer, parameter :: sim_specEpRatio=1.213


  !DEV : The Neitest config file in FLASH2 is hardcoded to use the single 
  !gamma EOS implementation.  We cannot do this in FLASH3 (and rightly so) 
  !because of a conflict between single gamma and multispecies.
  !For Ionize, we use multi gamma EOS implementation in FLASH3 and assign the same 
  !gamma value to each species.  If different gamma values are required for 
  !each species then the only file that needs to be modified is Simulation_initSpecies.
  real, parameter :: sim_singleGamma = 1.6667

  integer,dimension(SPEC_NUM) :: sim_specNumElect,sim_specElement
  logical,dimension(SPEC_NUM) :: sim_specSelected
  real, dimension(SPEC_NUM) :: sim_specAtomElem,sim_specAbundance
  character(len=2),dimension(SPEC_NUM) :: sim_specElementSymbol
  real,save, dimension(NSPECIES) :: sim_specXfrac

end module Simulation_speciesData
