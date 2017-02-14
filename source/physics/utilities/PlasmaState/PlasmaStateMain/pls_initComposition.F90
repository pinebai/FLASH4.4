!!****if* source/physics/utilities/PlasmaState/PlasmaStateMain/pls_initComposition
!!
!! NAME
!!
!!  pls_initComposition
!!
!! SYNOPSIS
!!
!!  call pls_initComposition()
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   No arguments
!!
!! AUTOGENROBODOC
!!
!!
!!***

subroutine pls_initComposition()
  use pls_interface, ONLY : pls_setSpeciesElementsData
  implicit none

  call pls_setSpeciesElementsData()
  return
end subroutine pls_initComposition
