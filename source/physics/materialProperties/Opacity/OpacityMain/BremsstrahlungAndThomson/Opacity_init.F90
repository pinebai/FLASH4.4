!!****if* source/physics/materialProperties/Opacity/OpacityMain/BremsstrahlungAndThomson/Opacity_init
!!
!! NAME
!!
!!  Opacity_init
!!
!!
!! SYNOPSIS
!!
!!  call Opacity_init()
!!
!! DESCRIPTION
!!
!! Initialiazed data for the Constcm2g opacity model using run time
!! parameters.
!!
!! ARGUMENTS
!!
!!  none
!!
!! PARAMETERS
!!
!!  useOpacity
!!  op_absorbScale
!!  op_emitScale
!!  op_transScale
!!
!!***
subroutine Opacity_init()
  use Opacity_data, ONLY : op_emitScale, op_transScale, &
       op_absorbScale, op_useOpacity
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none

  call RuntimeParameters_get ("useOpacity",   op_useOpacity)

  call RuntimeParameters_get("op_emitScale", op_emitScale)
  call RuntimeParameters_get("op_transScale", op_transScale)
  call RuntimeParameters_get("op_absorbScale", op_absorbScale)

end subroutine Opacity_init
