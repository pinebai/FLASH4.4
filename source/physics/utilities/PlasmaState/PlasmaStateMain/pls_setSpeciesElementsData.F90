!!****if* source/physics/utilities/PlasmaState/PlasmaStateMain/pls_setSpeciesElementsData
!!
!! NAME
!!
!!  pls_setSpeciesElementsData
!!
!! SYNOPSIS
!!
!!  call pls_setSpeciesElementsData ()
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
subroutine pls_setSpeciesElementsData ()

  use PlasmaState_data,                ONLY : pls_maxAtomicNumber,          &
                                          pls_maxNelementsPerSpecies,   &
                                          pls_maxNSpeciesPerElement,    &
                                          pls_totalSpecies,             &
                                          pls_totalElements,            &
                                          pls_elementNumberofSpecies,   &
                                          pls_element2AtomicNumber,     &
                                          pls_element2AtomicWeight,     &
                                          pls_elements2Species,         &
                                          pls_speciesNumberofElements,  &
                                          pls_speciesWeights,           &
                                          pls_species2Elements,         &
                                          pls_species2FractionElements

  use Multispecies_interface,      ONLY : Multispecies_getProperty
  use Driver_interface,            ONLY : Driver_abortFlash

  implicit none

# include "Flash.h"
# include "Multispecies.h"

  real, parameter :: zero  = 0.0

  integer :: element
  integer :: n,m
  integer :: nElements
  integer :: species
  integer :: Zval

  real    :: Aspecies

  integer, allocatable :: nAtomsUsed      (:)
  integer, allocatable :: nElementsUsed   (:)
  real,    allocatable :: numberFractions (:)
  real,    allocatable :: Zvalues         (:)
  real,    allocatable :: Avalues         (:)
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
      pls_totalSpecies = NSPECIES
  else
      call Driver_abortFlash ('[pls_initComposition] ERROR: no species found (NSPECIES =< 0)')
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
  allocate (avalues         (1:MS_MAXELEMS))
  allocate (numberFractions (1:MS_MAXELEMS))
  allocate (nAtomsUsed      (1:pls_maxAtomicNumber))

  nAtomsUsed                = 0
  pls_maxNelementsPerSpecies = 0

  do species = 1,pls_totalSpecies

     call Multispecies_getProperty (SPECIES_BEGIN - 1 + species , MS_NUMELEMS , nElements)

     if (pls_maxAtomicNumber < nElements) then
         call Driver_abortFlash ('[pls_initComposition] ERROR: Too many atomic elements / species')
     end if

     Zvalues = zero
     call Multispecies_getProperty (SPECIES_BEGIN - 1 + species , MS_ZELEMS   , Zvalues)

     do element = 1,nElements
        Zval = nint (Zvalues (element))
        if (Zval < 1 .or. Zval > pls_maxAtomicNumber) then
            call Driver_abortFlash ('[pls_initComposition] ERROR: Bad atomic number Z value')
        end if
        nAtomsUsed (Zval) = nAtomsUsed (Zval) + 1
     end do

     pls_maxNelementsPerSpecies = max (nElements , pls_maxNelementsPerSpecies)

  end do

  pls_totalElements         = count  (nAtomsUsed /= 0)
  pls_maxNspeciesPerElement = maxval (nAtomsUsed)
!
!
!    ...Allocate all arrays related to the species <-> element info.
!
!
  allocate (pls_speciesNumberofElements  (1:pls_totalSpecies))
  allocate (pls_speciesWeights           (1:pls_totalSpecies))
  allocate (pls_elementNumberofSpecies   (1:pls_totalElements))
  allocate (pls_element2AtomicNumber     (1:pls_totalElements))
  allocate (pls_element2AtomicWeight     (1:pls_totalElements)) !KW

  allocate (pls_elements2Species         (1:pls_totalElements, 1:pls_maxNspeciesPerElement))
  allocate (pls_species2Elements         (1:pls_totalSpecies , 1:pls_maxNelementsPerSpecies))
  allocate (pls_species2FractionElements (1:pls_totalSpecies , 1:pls_totalElements))
!
!
!    ...Inititalize the arrays with zeros.
!
!
  pls_elementNumberofSpecies     = 0
  pls_element2AtomicNumber       = 0
  pls_element2AtomicWeight       = -1 !KW
  pls_speciesNumberofElements    = 0
  pls_speciesWeights             = zero
  pls_elements2Species           = 0
  pls_species2Elements           = 0
  pls_species2FractionElements   = zero
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

  do Zval = 1,pls_maxAtomicNumber
     if (nAtomsUsed (Zval) /= 0) then
         element = element + 1
         pls_element2AtomicNumber   (element) = Zval
         pls_elementNumberofSpecies (element) = nAtomsUsed (Zval)
         nAtomsUsed (Zval) = element
     end if
  end do

  allocate (nElementsUsed (1:pls_totalElements))

  nElementsUsed = 0

  do species = 1,pls_totalSpecies

     call Multispecies_getProperty (SPECIES_BEGIN - 1 + species , A ,            Aspecies)
     call Multispecies_getProperty (SPECIES_BEGIN - 1 + species , MS_NUMELEMS  , nElements)
     call Multispecies_getProperty (SPECIES_BEGIN - 1 + species , MS_ZELEMS    , Zvalues)
     call Multispecies_getProperty (SPECIES_BEGIN - 1 + species , MS_AELEMS    , avalues)
     call Multispecies_getProperty (SPECIES_BEGIN - 1 + species , MS_FRACTIONS , numberFractions)

     pls_speciesNumberofElements (species) = nElements
     pls_speciesWeights          (species) = Aspecies

     do n = 1,nElements
        Zval = nint (Zvalues (n))
        element = nAtomsUsed (Zval)
        m = nElementsUsed (element) + 1
        nElementsUsed (element) = m
        pls_species2FractionElements (species,element) = numberFractions (n)
        pls_species2elements         (species,      n) = element
        pls_elements2Species         (element,      m) = species
!!$        pls_speciesElement2Avalues   (species,element) = avalues (n) !KW ??
        pls_element2AtomicWeight     (element)         = avalues (n)
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
end subroutine pls_setSpeciesElementsData
