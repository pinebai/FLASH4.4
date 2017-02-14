!!****if* source/physics/materialProperties/Opacity/OpacityMain/BremsstrahlungAndThomson/Opacity
!!
!! NAME
!!  Opacity
!!
!! SYNOPSIS
!!  call Opacity(real(in)  :: soln,
!!               real(in)  :: ngrp,
!!               real(out) :: opacityAbsorption,
!!               real(out) :: opacityEmission,
!!               real(out) :: opacityTransport )
!!
!! DESCRIPTION
!!
!!  Computes absorption, emission and transport opacities for a
!!  particular zone. 
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
       op_emitScale,  &
       op_transScale, &
       op_absorbScale
  implicit none
 
  real, intent(in), dimension(:) :: soln
  integer, intent(in) :: ngrp
  real, intent(out) :: opacityAbsorption
  real, intent(out) :: opacityEmission
  real, intent(out) :: opacityTransport

  real :: density, temperature, hydrogen, freeFreeConst
!
!
!   ...Check, if opacities need to be evaluated at all.
!
!

  freeFreeConst = 0.64e+23

  if (.not.op_useOpacity) then
     opacityAbsorption = 0.0
     opacityEmission   = 0.0
     opacityTransport  = 0.0
     return
  end if
!

  density = soln(DENS_VAR)
#ifdef TELE_VAR
  temperature = max(soln(TELE_VAR),10000.)
#else
  temperature = max(soln(TEMP_VAR),10000.)
#endif
  hydrogen = soln(H1_SPEC)

  opacityAbsorption =  op_absorbScale * freeFreeConst * density * (temperature)**(-3.5) * density 
  opacityEmission   =  op_emitScale * freeFreeConst * density * (temperature)**(-3.5) * density
  opacityTransport  =  op_transScale * (0.2 * (1.0 +hydrogen)) * density ! + freeFreeConst * density * (temperature)**(-3.5)) * density
  return
end subroutine Opacity

