!!****f* source/Simulation/Simulation_initSpecies
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


implicit none
end subroutine Simulation_initSpecies
