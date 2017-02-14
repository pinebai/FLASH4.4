!!****if* source/physics/materialProperties/Opacity/OpacityMain/BremsstrahlungAndThomson/Opacity_data
!!
!! NAME
!!
!!  Opacity_data
!!
!! SYNOPSIS
!!  use Opacity_data
!!
!! DESCRIPTION
!!
!!  Defines and stores local data for BremsstrahlungAndThomson implementation of the
!!  Opacity unit.
!!  
!!***
module Opacity_data
  
  implicit none

  logical, save :: op_useOpacity

  real, save :: op_emitScale
  real, save :: op_transScale
  real, save :: op_absorbScale

end module Opacity_data
