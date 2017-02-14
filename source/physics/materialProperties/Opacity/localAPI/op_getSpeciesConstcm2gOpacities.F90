!!****if* source/physics/materialProperties/Opacity/localAPI/op_getSpeciesConstcm2gOpacities
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

  implicit none

  integer, intent (in) :: species

  return
end subroutine op_getSpeciesConstcm2gOpacities
