!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Constcm2g/op_getSpeciesConstcm2gOpacities
!!
!! NAME
!!
!!  op_getSpeciesConstcm2gOpacities
!!
!! SYNOPSIS
!!
!!  call op_getSpeciesConstcm2gOpacities (integer (in) :: species,
!!                                        real    (in) :: dens)
!!
!! DESCRIPTION
!!
!!  Gets the constant absorption, emission and/or transport opacities for a
!!  particular species.
!!
!! ARGUMENTS
!!
!!   species : The species index
!!
!!***
subroutine op_getSpeciesConstcm2gOpacities (species)

  use Opacity_data,    ONLY : op_absorptionKind,     &
                              op_emissionKind,       &
                              op_transportKind,      &
                              op_speciesOpacities

  use op_constcm2gData, ONLY : op_absorptionConstcm2g, &
                               op_emissionConstcm2g,   &
                               op_transportConstcm2g

  implicit none
  
#include "Opacity.h"
#include "constants.h"
  
  integer, intent (in) :: species
!
!
!   ...Copy the stored values.
!
!
  if (op_absorptionKind (species) == OP_CONSTCM2G) then
      op_speciesOpacities (ABSORPTION,species) = op_absorptionConstcm2g (species)
  end if

  if (op_emissionKind (species) == OP_CONSTCM2G) then
      op_speciesOpacities (  EMISSION,species) = op_emissionConstcm2g   (species)
  end if

  if (op_transportKind (species) == OP_CONSTCM2G) then
      op_speciesOpacities ( TRANSPORT,species) = op_transportConstcm2g  (species)
  end if
!
!
!   ...Ready! 
!
!
  return
end subroutine op_getSpeciesConstcm2gOpacities
