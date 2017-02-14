!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Constant/op_constantData
!!
!! NAME
!!
!!  op_constantData
!!
!! SYNOPSIS
!!
!!  use op_constantData
!!
!! DESCRIPTION
!!
!!  Defines and stores the local data for the constant opacities.
!!  
!!***
module op_constantData
  
  implicit none

  logical, save :: op_initializedConstant = .false.

  real, allocatable, save :: op_absorptionConstant (:)
  real, allocatable, save :: op_emissionConstant   (:)
  real, allocatable, save :: op_transportConstant  (:)

end module op_constantData

