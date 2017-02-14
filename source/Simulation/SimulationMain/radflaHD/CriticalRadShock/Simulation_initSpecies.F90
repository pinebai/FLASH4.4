!!****if* source/Simulation/SimulationMain/radflaHD/CriticalRadShock/Simulation_initSpecies
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
use Multispecies_interface
use RuntimeParameters_interface, ONLY : RuntimeParameters_get
use Simulation_data, ONLY: sim_gamma
implicit none
#include "Multispecies.h"
#include "Flash.h"
#include "Eos.h"

call RuntimeParameters_get("gamma",sim_gamma) ! use for "ion" part of matter Eos 

call Multispecies_setProperty(H1_SPEC,A,1.)
call Multispecies_setProperty(H1_SPEC,Z,1.)
call Multispecies_setProperty(H1_SPEC,EB,0.)
call Multispecies_setProperty(H1_SPEC,MS_EOSTYPE,EOS_GAM)
call Multispecies_setProperty(H1_SPEC, GAMMA, sim_gamma)

!!$call Multispecies_setProperty(HE4_SPEC,A,4.)
!!$call Multispecies_setProperty(HE4_SPEC,Z,2.)
!!$call Multispecies_setProperty(HE4_SPEC,EB,28.29603)
!!$call Multispecies_setProperty(HE4_SPEC,MS_EOSTYPE,EOS_GAM)
!!$call Multispecies_setProperty(HE4_SPEC, GAMMA,  1.66666666667)
!!$
!!$call Multispecies_setProperty(C12_SPEC,A,12.)
!!$call Multispecies_setProperty(C12_SPEC,Z,6.)
!!$call Multispecies_setProperty(C12_SPEC,EB,92.16294)
!!$call Multispecies_setProperty(C12_SPEC,MS_EOSTYPE,EOS_GAM)
!!$call Multispecies_setProperty(C12_SPEC, GAMMA,  1.66666666667)
!!$
!!$call Multispecies_setProperty(N14_SPEC,A,14.)
!!$call Multispecies_setProperty(N14_SPEC,Z,7.)
!!$call Multispecies_setProperty(N14_SPEC,EB,104.65998)
!!$call Multispecies_setProperty(N14_SPEC,MS_EOSTYPE,EOS_GAM)
!!$call Multispecies_setProperty(N14_SPEC, GAMMA,  1.66666666667)
!!$
!!$call Multispecies_setProperty(O16_SPEC,A,16.)
!!$call Multispecies_setProperty(O16_SPEC,Z,8.)
!!$call Multispecies_setProperty(O16_SPEC,EB,127.62093)
!!$call Multispecies_setProperty(O16_SPEC,MS_EOSTYPE,EOS_GAM)
!!$call Multispecies_setProperty(O16_SPEC, GAMMA,  1.66666666667)

end subroutine Simulation_initSpecies
