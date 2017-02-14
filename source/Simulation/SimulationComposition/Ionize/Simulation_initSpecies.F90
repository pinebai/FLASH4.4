!!****if* source/Simulation/SimulationComposition/Ionize/Simulation_initSpecies
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
!!  This routine will initialize the species and species values needed
!!  for setups that use the source term for Ionization, there is another similar routine 
!!  in place for Burn networks. The setups that want to use the multispecies
!!  capabilities of the code for something other than nuclear burning and ionization
!!  should include their own custom implementation of this routine
!!
!!
!!
!! ARGUMENTS
!!
!!  There are no arguments in this subroutine
!!
!!***

subroutine Simulation_initSpecies()

  use Simulation_speciesData, ONLY : SPEC_NUM,sim_specEpRatio, sim_specNumElect,&
       sim_specElement,sim_specSelected,sim_specAtomElem,sim_specAbundance,&
       sim_specElementSymbol, sim_singleGamma
  use Driver_interface, ONLY : Driver_abortFlash
  use Multispecies_interface, ONLY : Multispecies_setProperty
  use Simulation_interface, ONLY : Simulation_mapStrToInt
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"
  

  real, save :: emass,pmass,e_abar, ave_hmass,ave1,ave2,electrons,zbar
  integer,save :: isotope
  
  integer :: numSpec,i, j
  integer, parameter :: SPEC_UNIT=20
  real :: bindEnergy

  ave1=0.0; ave2=0.0

  call PhysicalConstants_get("proton mass",pmass)
  call PhysicalConstants_get("electron mass",emass)
  e_abar=emass/pmass

  open(unit=SPEC_UNIT,file="SpeciesList.txt")
  
  !! First check for the number of elements to be included 
  !! the array numElect reads in the number of electrons in the 
  !! atom of the element under consideration, and aNumElem gets the
  !! atomic number. 
  !! We expect the Flash.h file to have created numElect+1 placeholders
  !! in the Species list for each included element. We also expect all 
  !! these placeholders to be contiguous. So for example if "he" is included
  !! then its species index is expeted to be 3 (after hydrogen and electron)
  !! the fourth and fifth species indices should be occupied by the two 
  !! electrons associated with it. If the next and last 
  !! included element is oxygen then we expect the species list in Flash.h 
  !! to read as follows :
  !! #define H_SPEC 1
  !! #define E_SPEC 2         !In our Flash.h we will have E_SPEC 1, H_SPEC 2
  !! #define HE_SPEC 3
  !! #define HE1_SPEC 4
  !! #define HE2_SPEC 5
  !! #define O_SPEC 6
  !! #define O1_SPEC 7
  !! #define O2_SPEC 8
  !! #define O3_SPEC 9
  !! #define O4_SPEC 10
  !! #define O5_SPEC 11
  !! #define O6_SPEC 12
  !! #define O7_SPEC 13
  !! #define O8_SPEC 14
  !! This is because O has 8 electrons. NSPECIES here is 14, even though the
  !! number of elements included is 3 (hydrogen, helium and oxygen)
  !! All elements, their electrons, atomic number and abundacnes etc
  !! are stored in file "SpeciesList.txt"
  numSpec=2   !CD:  I think this is because Hydrogen and Electron are not in SpeciesList.txt.
  do i = 1,SPEC_NUM
     read(SPEC_UNIT,*)sim_specElementSymbol(i),sim_specNumElect(i),sim_specAtomElem(i),sim_specAbundance(i)

     !! DEV: Expect this call to return the number correponding 
     !!      to the element, not one of its electrons. From the example
     !!      above, for he it should return 3, for oxygen 6 and for all
     !!      other elements NONEXISTENT

     call Simulation_mapStrToInt(sim_specElementSymbol(i),sim_specElement(i),MAPBLOCK_UNK)

     print *, "Symbol: ", sim_specElementSymbol(i), " has position:", sim_specElement(i), "in unk"

     sim_specSelected(i)= (sim_specElement(i) /= NONEXISTENT)
     if(sim_specSelected(i)) then
        numSpec=numSpec+sim_specNumElect(i)+1
     else  !CD:  What is this else clause for???
        ave1=ave1+sim_specAbundance(i)
        ave2=ave2+sim_specAbundance(i)*sim_specAtomElem(i)
     end if
  end do

  close(unit=SPEC_UNIT)

  if(numSpec/=NSPECIES)call Driver_abortFlash("Simulation_initSpecies:number of species does not match specifications")

  ave_hmass = 2.0*(1.0+sim_specEpRatio*e_abar+ave2)/(1.0+sim_specEpRatio+ave1)

  !! Now that all the elements to be included have been found, add the 
  !! elements and their electrons into the multifluids database

  isotope=ELEC_SPEC
  zbar = 0.0
  electrons = 1.0
  bindEnergy=1.0
  call Multispecies_setProperty(isotope, A, e_abar)
  call Multispecies_setProperty(isotope, Z, zbar)
  call Multispecies_setProperty(isotope, E, electrons)
  call Multispecies_setProperty(isotope, EB, bindEnergy)
  call Multispecies_setProperty(isotope, GAMMA, sim_singleGamma)

  
  isotope=H_SPEC
  zbar=1.0
  electrons=0.0
  call Multispecies_setProperty(isotope, A, ave_hmass)
  call Multispecies_setProperty(isotope, Z, zbar)
  call Multispecies_setProperty(isotope, E, electrons)
  call Multispecies_setProperty(isotope, EB, bindEnergy)
  call Multispecies_setProperty(isotope, GAMMA, sim_singleGamma)
  
  do i = 1,SPEC_NUM
     if(sim_specSelected(i))then
        isotope=sim_specElement(i)
        do j=0,sim_specNumElect(i)
           zbar=sim_specNumElect(i)*1.0
           electrons=zbar-j*1.0
           call Multispecies_setProperty(isotope, A, sim_specAtomElem(i))
           call Multispecies_setProperty(isotope, Z, zbar)
           call Multispecies_setProperty(isotope, E, electrons)
           call Multispecies_setProperty(isotope, EB, bindEnergy)
           call Multispecies_setProperty(isotope, GAMMA, sim_singleGamma)
           isotope=isotope+1
        end do
     end if
  end do
  
  
  !Calculate sim_specXfrac() for each ion.
  !sim_specXfrac = mass fraction / population fraction
  call sim_initFrac()

end subroutine Simulation_initSpecies
