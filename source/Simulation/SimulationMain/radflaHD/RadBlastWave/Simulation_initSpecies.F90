!!****if* source/Simulation/SimulationMain/radflaHD/RadBlastWave/Simulation_initSpecies
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
implicit none
#include "Multispecies.h"
#include "Flash.h"
#include "Eos.h"

call Multispecies_setProperty(H1_SPEC,A,2.0)
call Multispecies_setProperty(H1_SPEC,Z,1.0)
call Multispecies_setProperty(H1_SPEC,EB,0.)
call Multispecies_setProperty(H1_SPEC,MS_EOSTYPE,EOS_GAM)
call Multispecies_setProperty(H1_SPEC, GAMMA, 1.66666666667)

!!$call Multispecies_setProperty(HE3_SPEC,A,3.)
!!$call Multispecies_setProperty(HE3_SPEC,Z,2.)
!!$call Multispecies_setProperty(HE3_SPEC,EB,7.71819)
!!$call Multispecies_setProperty(HE3_SPEC,MS_EOSTYPE,EOS_GAM)
!!$call Multispecies_setProperty(HE3_SPEC, GAMMA, 1.66666666667)
!!$
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
!!$
!!$call Multispecies_setProperty(NE20_SPEC,A,20.)
!!$call Multispecies_setProperty(NE20_SPEC,Z,10.)
!!$call Multispecies_setProperty(NE20_SPEC,EB,160.64788)
!!$call Multispecies_setProperty(NE20_SPEC,MS_EOSTYPE,EOS_GAM)
!!$call Multispecies_setProperty(NE20_SPEC, GAMMA,  1.66666666667)
!!$
!!$call Multispecies_setProperty(MG24_SPEC,A,24.)
!!$call Multispecies_setProperty(MG24_SPEC,Z,12.)
!!$call Multispecies_setProperty(MG24_SPEC,EB,198.25790)
!!$call Multispecies_setProperty(MG24_SPEC,MS_EOSTYPE,EOS_GAM)
!!$call Multispecies_setProperty(MG24_SPEC, GAMMA,  1.66666666667)
!!$
!!$call Multispecies_setProperty(SI28_SPEC,A,28.)
!!$call Multispecies_setProperty(SI28_SPEC,Z,14.)
!!$call Multispecies_setProperty(SI28_SPEC,EB,236.53790)
!!$call Multispecies_setProperty(SI28_SPEC,MS_EOSTYPE,EOS_GAM)
!!$call Multispecies_setProperty(SI28_SPEC, GAMMA,  1.66666666667)
!!$
!!$call Multispecies_setProperty(S32_SPEC,A,32.)
!!$call Multispecies_setProperty(S32_SPEC,Z,16.)
!!$call Multispecies_setProperty(S32_SPEC,EB,271.78250)
!!$call Multispecies_setProperty(S32_SPEC,MS_EOSTYPE,EOS_GAM)
!!$call Multispecies_setProperty(S32_SPEC, GAMMA,  1.66666666667)
!!$
!!$call Multispecies_setProperty(AR36_SPEC,A,36.)
!!$call Multispecies_setProperty(AR36_SPEC,Z,18.)
!!$call Multispecies_setProperty(AR36_SPEC,EB,306.72020)
!!$call Multispecies_setProperty(AR36_SPEC,MS_EOSTYPE,EOS_GAM)
!!$call Multispecies_setProperty(AR36_SPEC, GAMMA,  1.66666666667)
!!$
!!$call Multispecies_setProperty(CA40_SPEC,A,40.)
!!$call Multispecies_setProperty(CA40_SPEC,Z,20.)
!!$call Multispecies_setProperty(CA40_SPEC,EB,342.05680)
!!$call Multispecies_setProperty(CA40_SPEC,MS_EOSTYPE,EOS_GAM)
!!$call Multispecies_setProperty(CA40_SPEC, GAMMA,  1.66666666667)
!!$
!!$call Multispecies_setProperty(TI44_SPEC,A,44.)
!!$call Multispecies_setProperty(TI44_SPEC,Z,22.)
!!$call Multispecies_setProperty(TI44_SPEC,EB,375.47720)
!!$call Multispecies_setProperty(TI44_SPEC,MS_EOSTYPE,EOS_GAM)
!!$call Multispecies_setProperty(TI44_SPEC, GAMMA,  1.66666666667)
!!$
!!$call Multispecies_setProperty(CR48_SPEC,A,48.)
!!$call Multispecies_setProperty(CR48_SPEC,Z,24.)
!!$call Multispecies_setProperty(CR48_SPEC,EB,411.46900)
!!$call Multispecies_setProperty(CR48_SPEC,MS_EOSTYPE,EOS_GAM)
!!$call Multispecies_setProperty(CR48_SPEC, GAMMA,  1.66666666667)
!!$
!!$call Multispecies_setProperty(FE52_SPEC,A,52.)
!!$call Multispecies_setProperty(FE52_SPEC,Z,26.)
!!$call Multispecies_setProperty(FE52_SPEC,EB,447.70800)
!!$call Multispecies_setProperty(FE52_SPEC,MS_EOSTYPE,EOS_GAM)
!!$call Multispecies_setProperty(FE52_SPEC, GAMMA,  1.66666666667)
!!$
!!$call Multispecies_setProperty(FE54_SPEC,A,54)
!!$call Multispecies_setProperty(FE54_SPEC,Z,26)
!!$call Multispecies_setProperty(FE54_SPEC,EB,471.76300)
!!$call Multispecies_setProperty(FE54_SPEC,MS_EOSTYPE,EOS_GAM)
!!$call Multispecies_setProperty(FE54_SPEC, GAMMA,  1.66666666667)
!!$
!!$call Multispecies_setProperty(NI56_SPEC,A,56.)
!!$call Multispecies_setProperty(NI56_SPEC,Z,28.)
!!$call Multispecies_setProperty(NI56_SPEC,EB,484.00300)
!!$call Multispecies_setProperty(NI56_SPEC,MS_EOSTYPE,EOS_GAM)
!!$call Multispecies_setProperty(NI56_SPEC, GAMMA,  1.66666666667)
!!$
!!$call Multispecies_setProperty(PROT_SPEC,A,1.)
!!$call Multispecies_setProperty(PROT_SPEC,Z,1.)
!!$call Multispecies_setProperty(PROT_SPEC,EB,0.)
!!$call Multispecies_setProperty(PROT_SPEC,MS_EOSTYPE,EOS_GAM)
!!$call Multispecies_setProperty(PROT_SPEC, GAMMA,  1.66666666667)
!!$
!!$call Multispecies_setProperty(NEUT_SPEC,A,1.)
!!$call Multispecies_setProperty(NEUT_SPEC,Z,0.)
!!$call Multispecies_setProperty(NEUT_SPEC,EB,0.)
!!$call Multispecies_setProperty(NEUT_SPEC,MS_EOSTYPE,EOS_GAM)
!!$call Multispecies_setProperty(NEUT_SPEC, GAMMA,  1.66666666667)

end subroutine Simulation_initSpecies
