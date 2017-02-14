!!****if* source/physics/materialProperties/Opacity/OpacityMain/Constcm2g/Opacity_data
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
!!  Defines and stores local data for Constcm2g implementation of the
!!  Opacity unit.
!!  
!!***
module Opacity_data
  
  implicit none

  logical, save :: op_useOpacity

  real, save :: op_emitConst
  real, save :: op_transConst
  real, save :: op_absorbConst

end module Opacity_data
