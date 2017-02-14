!!****if* source/physics/materialProperties/Opacity/OpacityMain/Constcm2g/Opacity
!!
!! NAME
!!  Opacity
!!
!! SYNOPSIS
!!  call Opacity(real(in)  :: soln
!!               real(in)  :: ngrp
!!               real(out) :: opacityAbsorption,
!!               real(out) :: opacityEmission,
!!               real(out) :: opacityTransport )
!!
!! DESCRIPTION
!!
!!  Computes absorption, emission and transport opacities for a
!!  particular zone. In this case the opacities are all set to a
!!  constant value defined in the flash.par file.
!!
!! ARGUMENTS
!!
!!   soln : The solution vector for the cell
!!   ngrp : The group number
!!
!!   opacityAbsorption : the absorption opacity (in 1/cm)
!!
!!   opacityEmission   : the emission opacity (in 1/cm)
!!
!!   opacityTransport  : the transport opacity which is used to compute the
!!                       diffusion coefficient (in 1/cm)
!!
!!***

#include "Flash.h"

subroutine Opacity(soln, ngrp, opacityAbsorption, opacityEmission, opacityTransport)
  use Opacity_data, ONLY: op_useOpacity, &
       op_emitConst,  &
       op_transConst, &
       op_absorbConst
  implicit none
 
  real, intent(in), dimension(:) :: soln
  integer, intent(in) :: ngrp
  real, intent(out) :: opacityAbsorption
  real, intent(out) :: opacityEmission
  real, intent(out) :: opacityTransport

  real :: density
!
!
!   ...Check, if opacities need to be evaluated at all.
!
!
  if (.not.op_useOpacity) then
     opacityAbsorption = 0.0
     opacityEmission   = 0.0
     opacityTransport  = 0.0
     return
  end if
!

  density = soln(DENS_VAR)

  opacityAbsorption = op_absorbConst * density
  opacityEmission   = op_emitConst   * density
  opacityTransport  = op_transConst  * density
  return
end subroutine Opacity

