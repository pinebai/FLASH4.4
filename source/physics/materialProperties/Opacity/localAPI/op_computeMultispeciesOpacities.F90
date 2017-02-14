!!****if* source/physics/materialProperties/Opacity/localAPI/op_computeMultispeciesOpacities
!!
!! NAME
!!
!!  op_computeMultispeciesOpacities
!!
!! SYNOPSIS
!!
!!  call op_computeMultispeciesOpacities (real (out) :: opacityAbsorption,
!!                                        real (out) :: opacityEmission,
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

  implicit none

  real,    intent (out) :: opacityAbsorption
  real,    intent (out) :: opacityEmission
  real,    intent (out) :: opacityTransport

  opacityAbsorption = 0.0
  opacityEmission   = 0.0
  opacityTransport  = 0.0

  return
end subroutine op_computeMultispeciesOpacities
