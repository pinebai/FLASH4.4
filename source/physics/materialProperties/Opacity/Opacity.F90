!!****f* source/physics/materialProperties/Opacity/Opacity
!!
!! NAME
!!
!!  Opacity
!!
!! SYNOPSIS
!!
!!  call Opacity (real    (in)  :: soln(:),
!!                integer (in)  :: ngrp,
!!                real    (out) :: opacityAbsorption,
!!                real    (out) :: opacityEmission,
!!                real    (out) :: opacityTransport )
!!
!! DESCRIPTION
!!
!!  Computes absorption, emission and transport opacities for a particular cell.
!!  Stub implementation sets all opacities to zero.
!!
!! ARGUMENTS
!!
!!   soln              : The solution vector for the cell
!!   ngrp              : The energy group number
!!   opacityAbsorption : the absorption opacity (in 1/cm)
!!   opacityEmission   : the emission opacity (in 1/cm)
!!   opacityTransport  : the transport opacity (in 1/cm)
!!
!!***

subroutine Opacity (soln, ngrp, opacityAbsorption, opacityEmission, opacityTransport)

  implicit none
  
  integer, intent(in)  :: ngrp
  real,    intent(out) :: opacityAbsorption
  real,    intent(out) :: opacityEmission
  real,    intent(out) :: opacityTransport

  real,    intent(in), dimension (:) :: soln
!
!
!      ...Stub values.
!
!
  opacityAbsorption = 0.0
  opacityEmission   = 0.0
  opacityTransport  = 0.0

  return
end subroutine Opacity
