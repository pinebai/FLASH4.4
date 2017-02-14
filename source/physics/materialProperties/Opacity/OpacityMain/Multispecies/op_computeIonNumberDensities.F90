!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/op_computeIonNumberDensities
!!
!! NAME
!!
!!  op_computeIonNumberDensities
!!
!! SYNOPSIS
!!
!!  call op_computeIonNumberDensities ()
!!
!! DESCRIPTION
!!
!!  Computes the ion number densities for all the species in the current cell
!!  and the total sum over all species. The individual ion number densities for
!!  each species is given by the following expression:
!!
!!          # ions / cm^3  =  Mf * rho * Na / A
!!
!!  where Mf and A are the mass fraction and weight of the species and
!!  rho and Na are the total cell mass density and Avogadro's number.
!!
!! ARGUMENTS
!!
!!***
subroutine op_computeIonNumberDensities ()

  use Opacity_data,    ONLY : op_totalSpecies,             &
                              op_speciesWeights,           &
                              op_cellMassDensity,          &
                              op_cellSpeciesMassFractions, &
                              op_ionNumberDensity,         &
                              op_totalIonNumberDensity

  use op_numericsData, ONLY : zero,                        &
                              op_Avogadro

  implicit none
  
#include "Flash.h"  
#include "constants.h"

  real :: A
  real :: Mf
  real :: numberDensity
  real :: rho

  integer :: species
!
!
!   ...Loop over all species.
!
!
  rho = op_cellMassDensity                                                   ! in g/cm^3

  op_totalIonNumberDensity = zero

  do species = 1,op_totalSpecies
     A                             = op_speciesWeights           (species)   ! in g/mol
     Mf                            = op_cellSpeciesMassFractions (species)   ! no dimensions
     numberDensity                 = rho * op_Avogadro / A                   ! in # ions/cm^3
     op_ionNumberDensity (species) = Mf * numberDensity

     op_totalIonNumberDensity      = op_totalIonNumberDensity + op_ionNumberDensity (species)
  end do
!
!
!   ...Ready! 
!
!
  return
end subroutine op_computeIonNumberDensities
