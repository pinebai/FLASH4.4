!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/op_setSpeciesElementsData
!!
!! NAME
!!
!!  op_setSpeciesElementsData
!!
!! SYNOPSIS
!!
!!  call op_setSpeciesElementsData ()
!!
!! DESCRIPTION
!!
!!  This routine sets all the data related to the species and their constituting elements.
!!  Several arrays are allocated and filled with the corresponding species and elements
!!  data. These arrays contain a complete picture of each species in terms of elements
!!  and vice versa.
!!
!! ARGUMENTS
!!
!!***
subroutine op_setSpeciesElementsData ()

  use Opacity_data,                ONLY : op_maxAtomicNumber,          &
                                          op_maxNelementsPerSpecies,   &
                                          op_maxNSpeciesPerElement,    &
                                          op_totalSpecies,             &
                                          op_totalElements,            &
                                          op_elementNumberofSpecies,   &
                                          op_element2AtomicNumber,     &
                                          op_elements2Species,         &
                                          op_speciesNumberofElements,  &
                                          op_speciesWeights,           &
                                          op_speciesLowTempCutoff,     &
                                          op_species2Elements,         &
                                          op_species2FractionElements

  use op_numericsData,             ONLY : zero

  use Multispecies_interface,      ONLY : Multispecies_getProperty
  use Driver_interface,            ONLY : Driver_abortFlash

  implicit none

# include "Flash.h"
# include "Multispecies.h"

  integer :: element
  integer :: n,m
  integer :: nElements
  integer :: species
  integer :: Zval

  real    :: Aspecies
  real    :: lowTempCutoff

  integer, allocatable :: nAtomsUsed      (:)
  integer, allocatable :: nElementsUsed   (:)
  real,    allocatable :: numberFractions (:)
  real,    allocatable :: Zvalues         (:)
!
!
!    ...Set the total # of species.
!
!       The concept of 'species' is an abstract one and is defined for each species
!       through its atomic element composition (which atomic elements and number
!       fraction composition).
!
!       The connection between the species and the corresponding atomic elements is
!       currently handled by the Multispecies unit and is not fully mature yet.
!
!
  if (NSPECIES > 0) then
      op_totalSpecies = NSPECIES
  else
      call Driver_abortFlash ('[Opacity_init] ERROR: no species found (NSPECIES =< 0)')
  end if
!
!
!    ...Determine the total # of different atomic elements present in the collection
!       of all species. For each species its atomic element composition (which atomic
!       elements and number fraction composition) is retrieved from the Multispecies unit.
!
!       Meaning of intermediately allocated arrays at this stage:
!
!               Zvalues = will contain the atomic number collection for each species
!       numberFractions = will contain the number fractions of each atomic element for
!                         each species
!            nAtomsUsed = will indicate which and how many of each atomic element is
!                         present in the collection of all species.
!
!
  allocate (Zvalues         (1:MS_MAXELEMS))
  allocate (numberFractions (1:MS_MAXELEMS))
  allocate (nAtomsUsed      (1:op_maxAtomicNumber))

  nAtomsUsed                = 0
  op_maxNelementsPerSpecies = 0

  do species = 1,op_totalSpecies

     call Multispecies_getProperty (SPECIES_BEGIN - 1 + species , MS_NUMELEMS , nElements)

     if (op_maxAtomicNumber < nElements) then
         call Driver_abortFlash ('[Opacity_init] ERROR: Too many atomic elements / species')
     end if

     Zvalues = zero
     call Multispecies_getProperty (SPECIES_BEGIN - 1 + species , MS_ZELEMS   , Zvalues)

     do element = 1,nElements
        Zval = nint (Zvalues (element))
        if (Zval < 1 .or. Zval > op_maxAtomicNumber) then
            call Driver_abortFlash ('[Opacity_init] ERROR: Bad atomic number Z value')
        end if
        nAtomsUsed (Zval) = nAtomsUsed (Zval) + 1
     end do

     op_maxNelementsPerSpecies = max (nElements , op_maxNelementsPerSpecies)

  end do

  op_totalElements         = count  (nAtomsUsed /= 0)
  op_maxNspeciesPerElement = maxval (nAtomsUsed)
!
!
!    ...Allocate all arrays related to the species <-> element info.
!
!
  allocate (op_speciesNumberofElements  (1:op_totalSpecies))
  allocate (op_speciesWeights           (1:op_totalSpecies))
  allocate (op_speciesLowTempCutoff     (1:op_totalSpecies))
  allocate (op_elementNumberofSpecies   (1:op_totalElements))
  allocate (op_element2AtomicNumber     (1:op_totalElements))

  allocate (op_elements2Species         (1:op_totalElements, 1:op_maxNspeciesPerElement))
  allocate (op_species2Elements         (1:op_totalSpecies , 1:op_maxNelementsPerSpecies))
  allocate (op_species2FractionElements (1:op_totalSpecies , 1:op_totalElements))
!
!
!    ...Inititalize the arrays with zeros.
!
!
  op_elementNumberofSpecies     = 0
  op_element2AtomicNumber       = 0
  op_speciesNumberofElements    = 0
  op_speciesWeights             = zero
  op_speciesLowTempCutoff       = zero
  op_elements2Species           = 0
  op_species2Elements           = 0
  op_species2FractionElements   = zero
!
!
!    ...Get the individual species <-> elements data.
!
!       The species are characterized by an index integer in Flash.h that was calculated from
!       the SPECIES_BEGIN index value. This might change in the future!
!
!       Meaning of intermediately allocated arrays at this stage:
!
!               nElementsUsed = will count how many of each atomic element has
!                               been used as we move through all the species
!
!
  element = 0

  do Zval = 1,op_maxAtomicNumber
     if (nAtomsUsed (Zval) /= 0) then
         element = element + 1
         op_element2AtomicNumber   (element) = Zval
         op_elementNumberofSpecies (element) = nAtomsUsed (Zval)
         nAtomsUsed (Zval) = element
     end if
  end do

  allocate (nElementsUsed (1:op_totalElements))

  nElementsUsed = 0

  do species = 1,op_totalSpecies

     call Multispecies_getProperty (SPECIES_BEGIN - 1 + species , A ,            Aspecies)
     call Multispecies_getProperty (SPECIES_BEGIN - 1 + species , MS_NUMELEMS  , nElements)
     call Multispecies_getProperty (SPECIES_BEGIN - 1 + species , MS_ZELEMS    , Zvalues)
     call Multispecies_getProperty (SPECIES_BEGIN - 1 + species , MS_FRACTIONS , numberFractions)
     call Multispecies_getProperty (SPECIES_BEGIN - 1 + species , MS_OPLOWTEMP , lowTempCutoff)

     op_speciesNumberofElements (species) = nElements
     op_speciesWeights          (species) = Aspecies
     op_speciesLowTempCutoff    (species) = lowTempCutoff         ! in Kelvin

     do n = 1,nElements
        Zval = nint (Zvalues (n))
        element = nAtomsUsed (Zval)
        m = nElementsUsed (element) + 1
        nElementsUsed (element) = m
        op_species2FractionElements (species,element) = numberFractions (n)
        op_species2elements         (species,      n) = element
        op_elements2Species         (element,      m) = species
     end do

  end do
!
!
!    ...Deallocate the temporary arrays.
!
!
  deallocate (Zvalues)
  deallocate (numberFractions)
  deallocate (nAtomsUsed)
  deallocate (nElementsUsed)
!
!
!    ...Ready!
!
!
  return
end subroutine op_setSpeciesElementsData
