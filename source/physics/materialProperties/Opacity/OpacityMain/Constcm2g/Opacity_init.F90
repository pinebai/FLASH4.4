!!****if* source/physics/materialProperties/Opacity/OpacityMain/Constcm2g/Opacity_init
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
!!  Initializes data for the Constcm2g opacity model using run time
!!  parameters.
!!
!! ARGUMENTS
!!
!!  none
!!
!! PARAMETERS
!!
!!  useOpacity
!!  op_absorbConst
!!  op_emitConst
!!  op_transConst
!!
!!***
subroutine Opacity_init()
  use Opacity_data, ONLY : op_emitConst, op_transConst, &
       op_absorbConst, op_useOpacity
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none

  call RuntimeParameters_get ("useOpacity",   op_useOpacity)

  call RuntimeParameters_get("op_emitConst", op_emitConst)
  call RuntimeParameters_get("op_transConst", op_transConst)
  call RuntimeParameters_get("op_absorbConst", op_absorbConst)

end subroutine Opacity_init
