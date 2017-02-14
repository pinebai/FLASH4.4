!!****if* source/physics/sourceTerms/Flame/FlameEffects/fl_effInit
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
!!  This version is stub-liek and does not do much.
!!
!! ARGUMENTS
!!
!!   No arguments
!!
!! HISTORY
!!
!! Dean Townsley 2008
!!***

subroutine fl_effInit()

  use fl_effData, only : fl_effsmlrho
  use RuntimeParameters_interface, only : RuntimeParameters_get

  implicit none

  ! for now just use global small density, could change in future
  call RuntimeParameters_get("smlrho", fl_effsmlrho)


end subroutine fl_effInit

