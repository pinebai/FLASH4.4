!!****f* source/Multispecies/Multispecies_getProperty
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

  implicit none


  integer, intent(in)       :: name, property
  real, intent(out)         :: value

  value = 0.0

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

subroutine Multispecies_getIntegerProperty(name, property, value)
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  integer, intent(in)       :: name, property
  integer, intent(out)      :: value

  value = -1


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
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  integer, intent(in)       :: name, property
  character(len=*), intent(out)      :: value

  value = ' '
end subroutine Multispecies_getStringProperty


subroutine Multispecies_getRealArrProperty(name, property, value)
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  integer, intent(in)           :: name, property
  real, intent(out)             :: value(:)
  value = 0.0

end subroutine Multispecies_getRealArrProperty


subroutine Multispecies_getIntArrProperty(name, property, value)
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  integer, intent(in)           :: name, property
  integer, intent(out)          :: value(:)
  value = 0

end subroutine Multispecies_getIntArrProperty
