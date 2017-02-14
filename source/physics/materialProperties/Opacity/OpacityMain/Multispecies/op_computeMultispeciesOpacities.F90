!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/op_computeMultispeciesOpacities
!!
!! NAME
!!
!!  op_computeMultispeciesOpacities
!!
!! SYNOPSIS
!!
!!  call op_computeMultispeciesOpacities (real (out) :: opacityAbsorption
!!                                        real (out) :: opacityEmission
!!                                        real (out) :: opacityTransport)
!!
!! DESCRIPTION
!!
!!  Computes the multispecies mass density scaled absorption, emission and transport
!!  opacities for a particular cell.
!!
!! ARGUMENTS
!!
!!   opacityAbsorption : The multispecies mass density scaled absorption opacity (in 1/cm)
!!   opacityEmission   : The multispecies mass density scaled   emission opacity (in 1/cm)
!!   opacityTransport  : The multispecies mass density scaled  transport opacity (in 1/cm)
!!
!!***
subroutine op_computeMultispeciesOpacities (opacityAbsorption, &
                                            opacityEmission,   &
                                            opacityTransport   )

  use Opacity_data,    ONLY : op_totalSpecies,          &
                              op_totalIonNumberDensity, &
                              op_ionNumberDensity,      &
                              op_cellMassDensity,       &
                              op_speciesOpacities

  use op_numericsData, ONLY : zero

  implicit none
  
#include "Flash.h"  
#include "Opacity.h"

  real, intent (out) :: opacityAbsorption
  real, intent (out) :: opacityEmission
  real, intent (out) :: opacityTransport

  integer :: species

  real    :: ionNumberDensity
  real    :: sumAbsorption
  real    :: sumEmission
  real    :: sumTransport
!
!
!   ...Loop over all species.
!
!
  sumAbsorption = zero
  sumEmission   = zero
  sumTransport  = zero

  do species = 1,op_totalSpecies

     ionNumberDensity = op_ionNumberDensity (species)

     sumAbsorption = sumAbsorption + ionNumberDensity * op_speciesOpacities (ABSORPTION,species)
     sumEmission   = sumEmission   + ionNumberDensity * op_speciesOpacities (  EMISSION,species)
     sumTransport  = sumTransport  + ionNumberDensity * op_speciesOpacities ( TRANSPORT,species)

  end do

  opacityAbsorption = op_cellMassDensity * (sumAbsorption / op_totalIonNumberDensity)
  opacityEmission   = op_cellMassDensity * (sumEmission   / op_totalIonNumberDensity)
  opacityTransport  = op_cellMassDensity * (sumTransport  / op_totalIonNumberDensity)
!
!
!   ...Ready! 
!
!
  return
end subroutine op_computeMultispeciesOpacities
