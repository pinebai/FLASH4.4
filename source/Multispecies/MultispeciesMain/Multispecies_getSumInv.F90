!!****if* source/Multispecies/MultispeciesMain/Multispecies_getSumInv
!!
!! NAME
!!
!!  Multispecies_getSumInv
!!
!!
!! SYNOPSIS
!!  Multispecies_getSumInv(integer(in) :: property,
!!                           real(out) :: value,
!!                   real(in),optional :: weights(:),
!!                integer(in),optional :: speciesMask(NSPECIES))
!!
!! DESCRIPTION
!!
!!  Given a property name and an array of weights, compute the
!!  weighted sum of the inverse of a chosen property.  
!!
!!  SumInv = sum over species( weight / propertyValue)
!!
!!  Weights should be an array of length equal to the number of
!!  species (NSPECIES).  If not provided, weights are assumed equal to
!!  one.
!!
!!  The default, is to return the weighted sum for
!!  all species.  A user, however, by passing in the optional argument
!!  speciesMask can choose to include only species of interest.
!!
!!  speciesMask is an array of integers of size Number of Species in length
!!  The mask works inclusively, that is the user specifies species that should be
!!  _included_ in the calculation.
!!
!!
!!  For example: in a given simulation if there are 4 defined species
!!  NI56, CO61, MG23, FE51 and the user wants to get a weighted sum of just
!!  NI56 and FE51 then the user would pass in the speciesMask arguement where
!!  speciesMask(1) = NI56
!!  speciesMask(2) = FE51
!!  speciesMask(3) = UNDEFINED_INT
!!  speciesMask(4) = UNDEFINED_INT
!!
!!  Order does not matter, however the order of the weights array MUST be the
!!  same as the mask.
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
!!    property - name of property as defined in Multispecies.h
!!    value - returned the sum of the inverses of a given property from routine
!!    weights - an optional array of length 1 or NSPECIES holding the weights for 
!!          each species.  If not provided, weight=1 is assumed.  If a weights
!!          array of size 1 is provided, its value is applied to all species.
!!    speciesMask - an optional array of length NSPECIES specifying the wanted 
!!              species for the average calculation
!!
!!  NOTES
!!
!!  Species properties are normally set in Simulation_initSpecies.
!!  The simulation's Config file defines the number and name of species, as in
!!     SPECIES AIR
!!     SPECIES SF6
!!
!! EXAMPLE
!!  subroutine Eos()
!!      use Eos_data, ONLY: 
!!      implicit none
!!  #include "Multispecies.h"
!!      .... declarations
!!      real, dimension(NSPECIES)      :: weights
!!      integer, dimension(NSPECIES)   :: mask
!!      .... executable statements
!!      call Multispecies_getSumInv(A,value,weights,mask)
!!      .....
!!  end subroutine Eos
!!
!!***

subroutine Multispecies_getSumInv(property, value, weights, speciesMask)

  use Multispecies_data  !, ONLY:  ms_isZero
  use Driver_interface, ONLY : Driver_abortFlash
  use Multispecies_interface, ONLY : Multispecies_getProperty

  implicit none

#include "Flash.h"
#include "Multispecies.h"

  integer, intent(in)                   :: property
  real, intent(out)                     :: value
  real, intent(in), optional            :: weights(:)
  integer, intent(in), optional         :: speciesMask(NSPECIES)

  real, dimension(NSPECIES)             :: weightsFull
  real                                  :: propVal, invPropVal
  integer                               :: numWeights, i, iMask

  value = 0.0

!! size of array of weights, mask, should be taken care of by compiler, 
!!    but double-check here.  Also generate a weights array "weightsFull" in case
!!    optional argument not given

  if (present(weights)) then
      numWeights = size(weights,1)
      if (numWeights == NSPECIES) then
         weightsFull = weights
      else if (numWeights == 1) then
         weightsFull = weights(1)
      else
         call Driver_abortFlash("Multispecies_getSumInv: invalid weights array")
      endif
   else
      weightsFull = 1.0
  endif

  if (present(speciesMask)) then
      if (size(speciesMask,1) /= NSPECIES ) then
         call Driver_abortFlash("Multispecies_getSumInv: mask array wrong size")
      endif
  endif

  if (.not. present(speciesMask)) then
     do i=1, NSPECIES 
        call Multispecies_getProperty(SPECIES_BEGIN + (i-1), property, propVal)
        if (ms_isZero(propVal))                                                 &
     &    call Driver_abortFlash("Multispecies_getSumInv: divide by ZERO property")
        invPropVal = 1.e0 / propVal
        value = value + weightsFull(i) * invPropVal
     end do

  else  ! there is a speciesMask

     do i=1, NSPECIES
        iMask = speciesMask(i)
        if (iMask /= UNDEFINED_INT) then
           call Multispecies_getProperty(iMask, property, propVal)

           if (ms_isZero(propVal))                                              &
     &      call Driver_abortFlash("Multispecies_getSumInv: divide by ZERO property")
           invPropVal = 1.e0 / propVal       
           value = value + weightsFull(i) * invPropVal
        end if

     end do 

  end if


end subroutine Multispecies_getSumInv
