!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Constant/op_getSpeciesConstantOpacities
!!
!! NAME
!!
!!  op_getSpeciesConstantOpacities
!!
!! SYNOPSIS
!!
!!  call op_getSpeciesConstantOpacities (integer (in) :: species)
!!
!! DESCRIPTION
!!
!!  Gets the constant absorption, emission and/or transport opacities for a
!!  particular species.
!!
!! ARGUMENTS
!!
!!   species   : The species index
!!   dens      : The mass density
!!
!!***
subroutine op_getSpeciesConstantOpacities (species, dens)

  use Opacity_data,    ONLY : op_absorptionKind,     &
                              op_emissionKind,       &
                              op_transportKind,      &
                              op_speciesOpacities

  use op_constantData, ONLY : op_absorptionConstant, &
                              op_emissionConstant,   &
                              op_transportConstant

  implicit none
  
#include "Opacity.h"
#include "constants.h"
  
  integer, intent (in) :: species
  real, intent (in) :: dens
!
!
!   ...Copy the stored values.
!
!
  if (op_absorptionKind (species) == OP_CONSTANT) then
      op_speciesOpacities (ABSORPTION,species) = op_absorptionConstant (species) / dens
  end if

  if (op_emissionKind (species) == OP_CONSTANT) then
      op_speciesOpacities (  EMISSION,species) = op_emissionConstant   (species) / dens
  end if

  if (op_transportKind (species) == OP_CONSTANT) then
      op_speciesOpacities ( TRANSPORT,species) = op_transportConstant  (species) / dens
  end if
!
!
!   ...Ready! 
!
!
  return
end subroutine op_getSpeciesConstantOpacities
