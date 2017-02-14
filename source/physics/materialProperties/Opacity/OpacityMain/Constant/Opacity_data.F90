!!****if* source/physics/materialProperties/Opacity/OpacityMain/Constant/Opacity_data
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
!!  Defines and stores the local data for Unit Opacity.
!!  
!!***
module Opacity_data
  
  implicit none

  logical, save :: op_useOpacity

  real, save :: op_emitConst
  real, save :: op_transConst
  real, save :: op_absorbConst

end module Opacity_data
