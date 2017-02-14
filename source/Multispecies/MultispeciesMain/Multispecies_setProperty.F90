!!****if* source/Multispecies/MultispeciesMain/Multispecies_setProperty
!!
!! NAME
!!
!!  Multispecies_setProperty
!!
!! 
!! SYNOPSIS
!!  Multispecies_setProperty(integer(in) :: name,
!!                           integer(in) :: property,
!!                           real(in)    :: value)
!!
!! DESCRIPTION
!!  Sets the property of the species name to the specified value.
!!  The name is an integer because it is defined in Flash.h based on Config data.
!!  
!!  The property is an integer because it is defined in Multispecies.h.
!!  description    property(defined as integer in Multispecies.h)
!!  --------------------------------------------------------------
!!  numTotal        A         Total number of protons and neutrons in nucleus
!!  numPositive     Z         Atomic number; number of protons in nucleus
!!  numNeutral      N         Number of neutrons
!!  numNegative     E         Number of electrons
!!  bindingEnergy   EB        Binding energy
!!  adiabatic index GAMMA     Ratio of heat capacities: Cp / Cv
!!  zMin            MS_ZMIN   a lower bound for Z
!!  
!!  ARGUMENTS
!!    name - name of species defined in Flash.h, e.g., NI56_SPEC
!!    property - name of property define as an integer
!!    value - value to set property to
!!
!!  NOTES
!!
!!  Species properties are normally set in Simulation_initSpecies.
!!  The simulation's Config file defines the number and name of species, as in
!!     SPECIES AIR
!!     SPECIES SF6
!!
!!
!!
!!***  

subroutine Multispecies_setProperty(name, property, value)


  use Multispecies_data  !, ONLY : ms_Array
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "Multispecies.h"


  integer, intent(in)           :: name, property
  real, intent(in)              :: value

  integer                       :: msindex

  call ms_mapMSIndex(name, msindex)

  if (property == A) then
     ms_Array(msindex)%numTotal = value
  else if (property == Z) then
     ms_Array(msindex)%numPositive = value
  else if (property == N) then
     ms_Array(msindex)%numNeutral = value
  else if (property == E) then
     ms_Array(msindex)%numNegative = value
  else if (property == EB) then
     ms_Array(msindex)%bindingEnergy = value
  else if (property == GAMMA) then
     ms_Array(msindex)%adiabaticIndex = value
  else if (property == MS_ZMIN) then
     ms_Array(msindex)%zMin = value
  else if (property == MS_OPLOWTEMP) then
     ms_Array(msindex)%opacityLowTemp = value
  else if (property == MS_EOSTYPE) then
     ms_Array(msindex)%eosType = value
  else if (property == MS_EOSSUBTYPE) then
     ms_Array(msindex)%eosSubType = value
  else
     call Driver_abortFlash("Error: Species property not found")
  end if



end subroutine Multispecies_setProperty

!!****if* source/Multispecies/MultispeciesMain/Multispecies_setProperty
!!
!! NAME
!!
!!  Multispecies_setIntegerProperty
!!
!! 
!! SYNOPSIS
!!  call Multispecies_setIntegerProperty(integer(in) :: name,
!!                                       integer(in) :: property,
!!                                       integer(in) :: value)
!!
!! DESCRIPTION
!!  Sets the integer property of the species name to the specified value.
!!  The name is an integer because it is defined in Flash.h based on Config data.
!!  
!!  The property is an integer because it is defined in Multispecies.h.
!!  description    property(defined as integer in Multispecies.h)
!!  --------------------------------------------------------------
!!  EOS type        MS_EOSTYPE  Type of EOS for this species
!!  
!!  ARGUMENTS
!!    name - name of species defined in Flash.h, e.g., NI56_SPEC
!!    property - name of property define as an integer
!!    value - value to set property to
!!
!!  NOTES
!!
!!  Species properties are normally set in Simulation_initSpecies.
!!  The simulation's Config file defines the number and name of species, as in
!!     SPECIES AIR
!!     SPECIES SF6
!!
!!
!!
!!***  

subroutine Multispecies_setIntegerProperty(name, property, value)


  use Multispecies_data  !, ONLY : ms_Array
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "Multispecies.h"


  integer, intent(in)           :: name, property
  integer, intent(in)           :: value

  integer                       :: msindex

  call ms_mapMSIndex(name, msindex)

  if (property == A) then
     ms_Array(msindex)%numTotal = value
  else if (property == Z) then
     ms_Array(msindex)%numPositive = value
  else if (property == N) then
     ms_Array(msindex)%numNeutral = value
  else if (property == E) then
     ms_Array(msindex)%numNegative = value
  else if (property == EB) then
     ms_Array(msindex)%bindingEnergy = value
  else if (property == GAMMA) then
     ms_Array(msindex)%adiabaticIndex = value
  else if (property == MS_ZMIN) then
     ms_Array(msindex)%zMin = value
  else if (property == MS_EOSTYPE) then
     ms_Array(msindex)%eosType = value
  else if (property == MS_EOSSUBTYPE) then
     ms_Array(msindex)%eossubType = value
  else if (property == MS_NUMELEMS) then
     ms_Array(msindex)%numElems = value
  else
     call Driver_abortFlash("Error: Species property not found")
  end if



end subroutine Multispecies_setIntegerProperty


!!****if* source/Multispecies/MultispeciesMain/Multispecies_setProperty
!!
!! NAME
!!
!!  Multispecies_setStringProperty
!!
!! 
!! SYNOPSIS
!!  call Multispecies_setStringProperty(integer(in) :: name,
!!                                       integer(in) :: property,
!!                                       character(len=*)(in) :: value)
!!
!! DESCRIPTION
!!  Sets the string property of the species name to the specified value.
!!  The name is an integer because it is defined in Flash.h based on Config data.
!!  
!!  The property is an integer because it is defined in Multispecies.h.
!!  description    property(defined as integer in Multispecies.h)
!!  --------------------------------------------------------------
!!  EOS electron file    MS_EOSELEFILE  name of file containing electron EOS table
!!  EOS ion file         MS_EOIONEFILE  name of file containing ion EOS table

!!  
!!  ARGUMENTS
!!    name - name of species defined in Flash.h, e.g., NI56_SPEC
!!    property - name of property define as an integer
!!    value - value to set property to
!!
!!  NOTES
!!
!!  Species properties are normally set in Simulation_initSpecies.
!!  The simulation's Config file defines the number and name of species, as in
!!     SPECIES AIR
!!     SPECIES SF6
!!
!!
!!
!!***  

subroutine Multispecies_setStringProperty(name, property, value)


  use Multispecies_data  !, ONLY : ms_Array
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "Multispecies.h"


  integer, intent(in)           :: name, property
  character(len=*), intent(in)           :: value

  integer                       :: msindex

  call ms_mapMSIndex(name, msindex)

  if (property == MS_EOSZFREEFILE) then
     ms_Array(msindex)%eosZfreetableFile = value
  else if (property == MS_EOSENERFILE) then
     ms_Array(msindex)%eosEnertableFile = value
  else if (property == MS_EOSPRESFILE) then
     ms_Array(msindex)%eosPrestableFile = value
  else if (property == MS_EOSGROUPNAME) then
     ms_Array(msindex)%eosGroupName = value
  else if (property == MS_EOSIONFILE) then
     print*,'WARNING: MS_EOISIONFILE is not valid in Multispecies_setStringProperty ',&
          'calls any more, IGNORED, pleade adjust your Simulation_initSpecies.'
  else if (property == MS_EOSELEFILE) then
     print*,'WARNING: MS_EOISELEFILE is not valid in Multispecies_setStringProperty ',&
          'calls any more, IGNORED, pleade adjust your Simulatele_initSpecies.'
  else
     call Driver_abortFlash("Error: String species property not found")
  end if



end subroutine Multispecies_setStringProperty




subroutine Multispecies_setRealArrProperty(name, property, value)
  use Multispecies_data  !, ONLY : ms_Array
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "Multispecies.h"

  integer, intent(in)           :: name, property
  real, intent(in)              :: value(:)
  integer                       :: msindex

  call ms_mapMSIndex(name, msindex)

  if (property == MS_AELEMS) then
     ms_Array(msindex)%aElems = value
  elseif (property == MS_ZELEMS) then
     ms_Array(msindex)%zElems = nint (value)
  elseif (property == MS_FRACTIONS) then
     ms_Array(msindex)%fractions = value
  else
     call Driver_abortFlash("Error: Species property not found")
  end if

end subroutine Multispecies_setRealArrProperty


subroutine Multispecies_setIntArrProperty(name, property, value)
  use Multispecies_data  !, ONLY : ms_Array
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "Multispecies.h"

  integer, intent(in)           :: name, property
  integer, intent(in)           :: value(:)
  integer                       :: msindex

  call ms_mapMSIndex(name, msindex)

  if (property == MS_ZELEMS) then
     ms_Array(msindex)%zElems = value
  else
     call Driver_abortFlash("Error: Species property not found")
  end if

end subroutine Multispecies_setIntArrProperty
