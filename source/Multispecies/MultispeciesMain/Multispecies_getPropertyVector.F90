!!****if* source/Multispecies/MultispeciesMain/Multispecies_getPropertyVector
!!
!! NAME
!!
!!  Multispecies_getPropertyVector
!!
!! 
!! SYNOPSIS
!!
!!  xall Multispecies_getPropertyVector(integer(in) :: property,
!!                                      real(out)   :: values)
!!
!! DESCRIPTION
!!
!!  Returns the values of a property of the species name in value.
!!  The name is an integer because it is defined in Flash.h based on Config data.
!!  
!!  description    property(defined as integer in Multispecies.h)
!!  --------------------------------------------------------------
!!  numTotal        A         Total number of protons and neutrons in nucleus
!!  numPositive     Z         Atomic number; number of protons in nucleus
!!  numNeutral      N         Number of neutrons
!!  numNegative     E         Number of electrons
!!  bindingEnergy   EB        Binding energy
!!  adiabatic index GAMMA     Ratio of heat capacities: Cp / Cv
!!  
!!  
!!  ARGUMENTS
!!    name - name of species defined in Flash.h, e.g., NI56_SPEC
!!    property - name of property define as an integer
!!    value - value of the returned property
!!
!!  NOTES
!!
!!  Species properties are normally set in Simulation_initSpecies.
!!  The simulation's Config file defines the number and name of species as
!!     SPECIES AIR
!!     SPECIES SF6
!!
!!***  

subroutine Multispecies_getPropertyVector(property, values)


  use Multispecies_data !, ONLY : ms_Array
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "Multispecies.h"


  integer, intent(in)       :: property
!!$  real, intent(out)         :: values(NSPECIES)
  real, intent(out)         :: values(:)


  if (property == A) then
     values = ms_Array%numTotal
  else if (property == Z) then
     values = ms_Array(:)%numPositive
  else if (property == N) then
     values = ms_Array(:)%numNeutral
  else if (property == E) then
     values = ms_Array(:)%numNegative
  else if (property == EB) then
     values = ms_Array(:)%bindingEnergy
  else if (property == GAMMA) then
     values = ms_Array(:)%adiabaticIndex
  else if (property == MS_EOSTYPE) then
     values = real(ms_Array(:)%eosType)
  else if (property == MS_EOSSUBTYPE) then
     values = real(ms_Array(:)%eosSubtype)
  else
     values = -1.0  ! default, to avoid compiler warnings with abort
     call Driver_abortFlash("Error: Species property not found")
  end if

end subroutine Multispecies_getPropertyVector

