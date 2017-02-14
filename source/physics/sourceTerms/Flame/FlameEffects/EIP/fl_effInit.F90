!!****if* source/physics/sourceTerms/Flame/FlameEffects/EIP/fl_effInit
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
!!
!! Initialize 
!!
!! ARGUMENTS
!!
!!   No arguments
!!
!!
!!
!!***

subroutine fl_effInit()

  use fl_effData
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

  call RuntimeParameters_get("flame_deltae",fl_effDeltae)
  call RuntimeParameters_get("ye_unburned",fl_eff_ye_u)
  call RuntimeParameters_get("sumyi_unburned",fl_eff_sumy_u)
  call RuntimeParameters_get("ye_burned",fl_eff_ye_b)
  call RuntimeParameters_get("sumyi_burned",fl_eff_sumy_b)

  call RuntimeParameters_get("smlrho", fl_effsmlrho)

end subroutine fl_effInit

