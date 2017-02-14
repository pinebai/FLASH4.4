!!****f* source/Multispecies/Multispecies_getSumFrac
!!
!! NAME
!!
!!  Multispecies_getSumFrac
!!
!!
!! SYNOPSIS
!!  Multispecies_getSumFrac(integer(in) :: property,
!!                            real(out) :: value,
!!                    real(in),optional :: weights(:),
!!                 integer(in),optional :: speciesMask(NSPECIES))
!!
!! DESCRIPTION
!!
!!  Given a property name and an array of weights, compute the
!!  fraction of the weighted sum of the chosen property. The sum is
!!  weighted by the total weight of the species (A).  
!!      SumFrac = sum over species (weights * propertyValue / A)
!!
!!  Weights should be an array of length equal to the number of
!!  species (NSPECIES).  If not provided, weights are assumed equal to
!!  one.
!!
!!  The default is to return the weighted sum fraction for
!!  all species.  A user, however, by passing in the optional argument
!!  speciesMask can choose to include only species of interest.
!!
!!  speciesMask is an array of integers of size NSPECIES.
!!  The mask works inclusively, that is, the user specifies species that
!!  should be _included_ in the calculation.
!!
!!  For example: in a given simulation there are 4 defined species
!!  NI56, CO61, MG23, FE51.  They are defined in Flash.h as NI56_SPEC,
!!  CO61_SPEC, MG23_SPEC, and FE51_SPEC respectively.  If the user
!!  wants to get a weighted sum of just NI56 and FE51 then the user would
!!  pass in the speciesMask argument where
!!
!!  speciesMask(1) = NI56_SPEC
!!  speciesMask(2) = FE51_SPEC
!!  speciesMask(3) = UNDEFINED_INT
!!  speciesMask(4) = UNDEFINED_INT
!!
!!  Order does not matter, however the order of the weights array MUST
!!  be the same as the mask.  For example, if the user wants NI56
!!  weighted twice as heavily as FE51 in the previous example, he
!!  would set
!!
!!   weights(1) = 2.0
!!   weights(2) = 1.0
!!   weights(3) = UNDEFINED_REAL
!!   weights(4) = UNDEFINED_REAL
!!
!!  The property is an integer because it is defined in Multispecies.h
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
!!
!!    property - name of property, as an integer defined in Multispecies.h
!!    value    - returned weighted sum fraction of a given property
!!    weights  - an optional array of length 1 or NSPECIES holding the weights for each
!!                species.  If not provided, weights=1 is assumed.
!!    speciesMask - optional array of length NSPECIES specifying the species
!!                     to be included in the average.
!!                  If not provided, all species are used.
!!
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
!! EXAMPLE
!!  subroutine Eos()
!!      use Eos_data, ONLY: 
!!      implicit none
!!#include "Multispecies.h"
!!      .... declarations
!!      real, dimension(NSPECIES)      :: weights
!!      integer, dimension(NSPECIES)   :: mask
!!      .... executable statements
!!      call Multispecies_getSumFrac(A,value,weights,mask)
!!      .....
!!   end subroutine Eos
!!
!! SEE ALSO
!!
!!   Multispecies_getSum
!!***

subroutine Multispecies_getSumFrac(property, value, weights, speciesMask)
  
  implicit none
#include "Flash.h"
  
  integer, intent(in)           :: property
  real, intent(out)             :: value
  real, intent(in), optional    :: weights(:)
  integer, intent(in), optional :: speciesMask(NSPECIES)

  value = 0.0
  
end subroutine Multispecies_getSumFrac
