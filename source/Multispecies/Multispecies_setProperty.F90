!!****f* source/Multispecies/Multispecies_setProperty
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

  implicit none

  integer, intent(in)           :: name, property
  real, intent(in)              :: value

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
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  integer, intent(in)           :: name, property
  integer, intent(in)           :: value

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
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  integer, intent(in)           :: name, property
  character(len=*), intent(in)           :: value

end subroutine Multispecies_setStringProperty




subroutine Multispecies_setRealArrProperty(name, property, value)
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  integer, intent(in)           :: name, property
  real, intent(in)              :: value(:)

end subroutine Multispecies_setRealArrProperty



subroutine Multispecies_setIntArrProperty(name, property, value)
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  integer, intent(in)           :: name, property
  integer, intent(in)           :: value(:)

end subroutine Multispecies_setIntArrProperty
