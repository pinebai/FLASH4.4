!!****if* source/Multispecies/MultispeciesMain/Multispecies_getProperty
!!
!! NAME
!!
!!  Multispecies_getProperty
!!
!! 
!! SYNOPSIS
!!
!!  Multispecies_getProperty(integer(in) :: name,
!!                           integer(in) :: property,
!!                           real(out)   :: value)
!!
!! DESCRIPTION
!!
!!  Returns the value of a property of the species name in value.
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
!!  zMin            MS_ZMIN   a lower bound for Z
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

subroutine Multispecies_getProperty(name, property, value)


  use Multispecies_data !, ONLY : ms_Array
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "Multispecies.h"


  integer, intent(in)       :: name, property
  real, intent(out)         :: value

  integer                   :: msindex


  call ms_mapMSIndex(name, msindex)

  if (property == A) then
     value = ms_Array(msindex)%numTotal
  else if (property == Z) then
     value = ms_Array(msindex)%numPositive
  else if (property == N) then
     value = ms_Array(msindex)%numNeutral
  else if (property == E) then
     value = ms_Array(msindex)%numNegative
  else if (property == EB) then
     value = ms_Array(msindex)%bindingEnergy
  else if (property == GAMMA) then
     value = ms_Array(msindex)%adiabaticIndex
  else if (property == MS_ZMIN) then
     value = ms_Array(msindex)%zMin
  else if (property == MS_OPLOWTEMP) then
     value = ms_Array(msindex)%opacityLowTemp
  else if (property == MS_EOSTYPE) then
     value = real(ms_Array(msindex)%eosType)
  else if (property == MS_EOSSUBTYPE) then
     value = real(ms_Array(msindex)%eosSubtype)
  else
     value = -1.0  ! default, to avoid compiler warnings with abort
     call Driver_abortFlash("Error: Species property not found")
  end if

end subroutine Multispecies_getProperty

!!****if* source/Multispecies/MultispeciesMain/Multispecies_getIntegerProperty
!!
!! NAME
!!
!!  Multispecies_getIntegerProperty
!!
!! 
!! SYNOPSIS
!!
!!  call Multispecies_getIntegerProperty(integer(in)  :: name,
!!                                       integer(in)  :: property,
!!                                       integer(out) :: value)
!!
!! DESCRIPTION
!!
!!  Returns the value of an integer property of the species name in value.
!!  The name is an integer because it is defined in Flash.h based on Config data.
!!  
!!  description    integer property(defined as integer in Multispecies.h)
!!  --------------------------------------------------------------
!!  EOS type        MS_EOSTYPE     Type of EOS for this species
!!  EOS subtype     MS_EOSSUBTYPE  Subtype of EOS for this species
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

subroutine Multispecies_getIntegerProperty(name, property, value)


  use Multispecies_data, ONLY : ms_Array
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "Multispecies.h"


  integer, intent(in)       :: name, property
  integer, intent(out)      :: value

  integer                   :: msindex


  call ms_mapMSIndex(name, msindex)

  if (property == A) then
     value = nint(ms_Array(msindex)%numTotal)
  else if (property == Z) then
     value = nint(ms_Array(msindex)%numPositive)
  else if (property == N) then
     value = nint(ms_Array(msindex)%numNeutral)
  else if (property == E) then
     value = nint(ms_Array(msindex)%numNegative)
  else if (property == MS_ZMIN) then
     value = nint(ms_Array(msindex)%zMin)
  else if (property == MS_EOSTYPE) then
     value = ms_Array(msindex)%eosType
  else if (property == MS_EOSSUBTYPE) then
     value = ms_Array(msindex)%eosSubtype
  else if (property == MS_NUMELEMS) then
     value = ms_Array(msindex)%numElems
  else
     value = -1  ! default, to avoid compiler warnings with abort
     call Driver_abortFlash("Error: Integer species property not found")
  end if

end subroutine Multispecies_getIntegerProperty


!!****if* source/Multispecies/MultispeciesMain/Multispecies_getStringProperty
!!
!! NAME
!!
!!  Multispecies_getStringProperty
!!
!! 
!! SYNOPSIS
!!
!!  call Multispecies_getStringProperty(integer(in)  :: name,
!!                                       integer(in)  :: property,
!!                                       character(len=*)(out) :: value)
!!
!! DESCRIPTION
!!
!!  Returns the value of an string property of the species name in value.
!!  The name is an integer because it is defined in Flash.h based on Config data.
!!  
!!  description    integer property(defined as integer in Multispecies.h)
!!  --------------------------------------------------------------
!!  EOS type        MS_EOSTYPE  Type of EOS for this species
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

subroutine Multispecies_getStringProperty(name, property, value)


  use Multispecies_data, ONLY : ms_Array
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "Multispecies.h"


  integer, intent(in)       :: name, property
  character(len=*), intent(out)      :: value

  integer                   :: msindex


  call ms_mapMSIndex(name, msindex)

  if (property == MS_EOSZFREEFILE) then
     value = ms_Array(msindex)%eosZfreetableFile
  else if (property == MS_EOSENERFILE) then
     value = ms_Array(msindex)%eosEnertableFile
  else if (property == MS_EOSPRESFILE) then
     value = ms_Array(msindex)%eosPrestableFile
  else if (property == MS_EOSGROUPNAME) then
     value = ms_Array(msindex)%eosGroupName
  else
     value = ' '  ! default, to avoid compiler warnings with abort
     call Driver_abortFlash("Error: String species property not found")
  end if

end subroutine Multispecies_getStringProperty


subroutine Multispecies_getRealArrProperty(name, property, value)
  use Multispecies_data  !, ONLY : ms_Array
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "Multispecies.h"

  integer, intent(in)           :: name, property
  real, intent(out)             :: value(:)
  integer                       :: msindex

  call ms_mapMSIndex(name, msindex)

  if (property == MS_AELEMS) then
     value = ms_Array(msindex)%aElems
  elseif (property == MS_ZELEMS) then
     value = real (ms_Array(msindex)%zElems)
  elseif (property == MS_FRACTIONS) then
     value = ms_Array(msindex)%fractions
  else
     value = 0.0
     call Driver_abortFlash("Error: Species property not found")
  end if

end subroutine Multispecies_getRealArrProperty


subroutine Multispecies_getIntArrProperty(name, property, value)
  use Multispecies_data  !, ONLY : ms_Array
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "Multispecies.h"

  integer, intent(in)           :: name, property
  integer, intent(out)          :: value(:)
  integer                       :: msindex

  call ms_mapMSIndex(name, msindex)

  if (property == MS_ZELEMS) then
     value = ms_Array(msindex)%zElems
  else
     value = 0
     call Driver_abortFlash("Error: Species property not found")
  end if

end subroutine Multispecies_getIntArrProperty
