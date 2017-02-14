!!****if* source/physics/sourceTerms/Flame/FlameEffects/BurnParametric/fl_effInit
!!
!! NAME
!!
!!  fl_effInit
!!
!! SYNOPSIS
!!
!!  call fl_effInit()
!!
!! DESCRIPTION
!!  This gets some runtime parameters for FlameEffects implementations.
!!
!! ARGUMENTS
!!
!!   No arguments
!!
!! HISTORY
!!
!! Dean Townsley 2008
!! Klaus Weide   2013
!!***


subroutine fl_effInit()

  use fl_effData, only : fl_effsmlrho, fl_effEosTol
  use RuntimeParameters_interface, only : RuntimeParameters_get

  implicit none

  ! for now just use global small density, could change in future
  call RuntimeParameters_get("smlrho", fl_effsmlrho)
  ! for now just use eos_tolerance runtime parameter shared with the Eos unit.
  call RuntimeParameters_get("eos_tolerance", fl_effEosTol)


end subroutine fl_effInit

