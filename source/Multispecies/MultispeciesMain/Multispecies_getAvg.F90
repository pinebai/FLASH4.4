!!****if* source/Multispecies/MultispeciesMain/Multispecies_getAvg
!!
!! NAME
!!
!!  Multispecies_getAvg
!!
!! 
!! SYNOPSIS
!!
!!  Multispecies_getAvg(     integer(in) :: property,
!!                             real(out) :: value,
!!                     real(in),optional :: weights(:),
!!                  integer(in),optional :: speciesMask(NSPECIES))
!!
!! DESCRIPTION
!!
!!  Gets the weighted average value of a property for requested
!!  species.  
!!    Avg = sum over species( weight * propertyValue) / (number of species)
!!
!!  The default is to return the average for all species.  A user,
!!  however, by passing in the optional argument speciesMask can
!!  choose to include only species of interest.
!!
!!  speciesMask is an array of integers of size Number of Species (NSPECIES) 
!!      in length.  The speciesMask works inclusively, that is the user 
!!      specifies species that should be _included_ in the average.
!!
!!  For example: in a given simulation there are 4 defined species
!!  NI56, CO61, MG23, FE51.  They are defined in Flash.h as NI56_SPEC,
!!  CO61_SPEC, MG23_SPEC, and FE51_SPEC respectively.  If the user
!!  wants to get an average of just NI56 and FE51 then the user would
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
!!    property    - name of property, as an integer defined in Multispecies.h
!!    value       - returned average value from routine
!!    weights     - an optional array of length 1 or NSPECIES holding the weights for 
!!                  each species.  If not provided, weight=1 is assumed.  If a weights
!!                  array of size 1 is provided, its value is applied to all species.
!!    speciesMask - an optional array of length NSPECIES specifying the wanted 
!!                     species for the average calculation
!!                  If not provided, all species are used.
!!
!! NOTES
!!  
!!  Species properties are normally set in Simulation_initSpecies.
!!  The simulation's Config file defines the number and name of species, as in
!!     SPECIES AIR
!!     SPECIES SF6
!!
!! EXAMPLE
!!
!!  subroutine Eos()
!!      use Eos_data, ONLY: 
!!      implicit none
!!#include "Multispecies.h"
!!      .... declarations
!!      real, dimension(NSPECIES)      :: weights
!!      integer, dimension(NSPECIES)   :: mask
!!      .... executable statements
!!      call Multispecies_getAvg(A,value,weights,mask)
!!      .....
!!  end subroutine Eos
!!
!!
!!
!!
!!***  

subroutine Multispecies_getAvg(property, value, weights, speciesMask)

  use Multispecies_interface, ONLY : Multispecies_getProperty
  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface, ONLY: Logfile_stamp

  implicit none
  
#include "Flash.h"
#include "Multispecies.h"

  integer, intent(in)               :: property
  real, intent(out)                 :: value
  real, intent(in), optional        :: weights(:)
  integer, intent(in), optional     :: speciesMask(NSPECIES)

  real                              :: tempSum, propVal
  real, dimension(NSPECIES)         :: weightsFull
  integer                           :: i, numWeights, iMask, iCount
  logical, save                     :: warningNeeded = .true.

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
         call Driver_abortFlash("Multispecies_getAvg: invalid weights array")
      endif
   else
      weightsFull = 1.0
  endif

  if (present(speciesMask)) then
      if (size(speciesMask,1) /= NSPECIES ) then
         call Driver_abortFlash("Multispecies_getAvg: mask array wrong size")
      endif
  endif

  tempSum = 0
  if (.not. present(speciesMask)) then
     do i=1, NSPECIES 
        
        call Multispecies_getProperty(SPECIES_BEGIN + (i-1), property, propVal)
        tempSum = tempSum + weightsFull(i) * propVal
        
     end do
  ! Calls of this routine should not happen when NSPECIES == 0 - KW
  !  Let's use regular warning systems  - LBR
#if (NSPECIES > 0) 
        value = tempSum / NSPECIES
#else
        value = 0.0
        if (warningNeeded) then
           call Logfile_stamp(NSPECIES, &
             '[Multispecies_getAvg] Return value set to 0 since NSPECIES = ]')
           write(*,*) &
             '[Multispecies_getAvg] Return value set to 0 since NSPECIES = ',NSPECIES 
           warningNeeded = .false.
        end if
#endif

  else  ! there is a speciesMask
     iCount = 0
     do i=1, NSPECIES
        iMask = speciesMask(i)
        if (iMask /= UNDEFINED_INT) then
           iCount = iCount + 1
           call Multispecies_getProperty(iMask, property, propVal)
           tempSum = tempSum + weightsFull(i) * propVal
        end if
     end do

     value = tempSum / iCount
     
  end if

end subroutine Multispecies_getAvg
