!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Constcm2g/op_constcm2gData
!!
!! NAME
!!
!!  op_constcm2gData
!!
!! SYNOPSIS
!!
!!  use op_constcm2gData
!!
!! DESCRIPTION
!!
!!  Defines and stores the local data for the constant specific
!!  opacities (in units of cm^2/g)
!!  
!!***
module op_constcm2gData
  
  implicit none

  logical, save :: op_initializedConstcm2g = .false.

  real, allocatable, save :: op_absorptionConstcm2g (:)
  real, allocatable, save :: op_emissionConstcm2g   (:)
  real, allocatable, save :: op_transportConstcm2g  (:)

end module op_constcm2gData
