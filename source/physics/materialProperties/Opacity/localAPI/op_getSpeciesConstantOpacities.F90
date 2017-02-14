!!****if* source/physics/materialProperties/Opacity/localAPI/op_getSpeciesConstantOpacities
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

  implicit none

  integer, intent (in) :: species
  real,    intent (in) :: dens

  return
end subroutine op_getSpeciesConstantOpacities
