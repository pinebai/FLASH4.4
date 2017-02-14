!!****if* source/physics/utilities/PlasmaState/PlasmaStateMain/PlasmaState_getComposition
!!
!! NAME
!!
!!  PlasmaState_getComposition
!!
!! SYNOPSIS
!!
!!  call PlasmaState_getComposition(real(OUT), dimension(:)  :: vecz,
!!                                  real(OUT), dimension(:)  :: veca,
!!                                  real(OUT), dimension(:)  :: vecfrac,
!!                                  integer(OUT)  :: nfrac,
!!                                  real(IN), dimension(:)  :: solnstate)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   vecz :       returns list of element Z values as a one-dimensional array
!!
!!   veca :       returns list of element A values as a one-dimensional array
!!
!!   vecfrac :    returns list of number fractions values as a one-dimensional array.
!!                They should add up to 1.
!!
!!   nfrac :      Numer of valid elements in the arrays that are returned.
!!
!!   solnstate :  solution state in one cell
!!
!!
!!
!!***


#include "Flash.h"
#include "constants.h"

subroutine PlasmaState_getComposition(vecZ, vecA, vecFrac, nFrac, solnState)
  use Multispecies_interface, ONLY : Multispecies_getPropertyVector,Multispecies_getSumInv, &
                                     Multispecies_getProperty
  use Eos_interface, ONLY : Eos_getAbarZbar
  use PlasmaState_data, ONLY : pls_meshMe, &
                               pls_totalElements, &
                               pls_elementNumberofSpecies, &
                               pls_element2AtomicNumber,   &
                               pls_element2AtomicWeight,   &
                               pls_elements2Species,       &
                               pls_species2FractionElements
  implicit none

#include "Multispecies.h"

  real, intent(OUT), dimension(:) :: vecZ, vecA, vecFrac
  integer, intent(OUT)            :: nFrac
  real, intent(IN),  dimension(:) :: solnState

  integer :: nout
  integer :: ie,s,species
  real    :: fraction,speciesMoleFraction
  real    :: abarInv, aspecies

  nout = min(size(vecZ),size(vecA), size(vecFrac))

#ifdef FLASH_MULTISPECIES
  if (pls_totalElements > 0) then
     nFrac = pls_totalElements
     call checkSpace()
     call Multispecies_getSumInv(A, abarInv,solnState(SPECIES_BEGIN:SPECIES_END))
     do ie = 1,nFrac
        vecZ(ie) = pls_element2AtomicNumber   (ie) ! = Zval
        vecA(ie) = pls_element2AtomicWeight   (ie)
        vecFrac(ie) = 0.0
        do s = 1, pls_elementNumberofSpecies (ie) ! = nAtomsUsed (Zval)
           species  = pls_elements2Species         (ie      ,s)
           fraction = pls_species2FractionElements (species,ie)
           call Multispecies_getProperty (SPECIES_BEGIN-1+species , A , Aspecies)
           speciesMoleFraction = solnState(SPECIES_BEGIN-1+species) / (abarInv * Aspecies)
           vecFrac(ie) = vecFrac(ie) + fraction * speciesMoleFraction
        end do
     end do


  else
     nFrac = NSPECIES
     call checkSpace()
     call Multispecies_getPropertyVector(A,vecA(:nFrac))
     call Multispecies_getPropertyVector(Z,vecZ(:nFrac))
     vecFrac(1:NSPECIES) = solnState(SPECIES_BEGIN:SPECIES_END)
  end if
#else
  nFrac = 1
  call checkSpace()
  call Eos_getAbarZbar(solnState,abar=vecA(1),zbar=vecZ(1))
  vecFrac(1) = 1.0

#endif

contains
  subroutine checkSpace()
    if (nout < nFrac) then
99     format('PlasmaState_getComposition: not enough space to return data for',&
          I4,'element(s)!') 
       if (pls_meshMe == MASTER_PE) print 99,nFrac
       call Driver_abortFlash('PlasmaState_getComposition: not enough space to return composition data')
    end if
  end subroutine checkSpace
end subroutine PlasmaState_getComposition
