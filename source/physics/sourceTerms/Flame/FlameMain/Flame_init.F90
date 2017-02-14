!!****if* source/physics/sourceTerms/Flame/FlameMain/Flame_init
!!
!! NAME
!!
!!  Flame_init
!!
!! SYNOPSIS
!!
!!  call Flame_init()
!!
!! DESCRIPTION
!!
!! Dean Townsley 2008
!! Initialize
!!
!!
!! ARGUMENTS
!!
!!   No arguments
!!
!!
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine Flame_init()

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Grid_interface, ONLY : Grid_getMinCellSize
  use Flame_data
  use fl_fsInterface, only : fl_fsGcMask, fl_fsInit
  use fl_effInterface, only : fl_effInit

  implicit none

  real :: dx

  ! get runtime parameters
  call RuntimeParameters_get("useFlame", fl_useFlame)
  call RuntimeParameters_get("fl_epsilon_0", fl_epsilon_0)
  call RuntimeParameters_get("fl_epsilon_1", fl_epsilon_1)
  call RuntimeParameters_get("fl_kpp_fact", fl_kpp_fact)
  call RuntimeParameters_get("fl_b", fl_b)

  call RuntimeParameters_get("fl_initProfileAdjustWidth", fl_initProfileAdjustWidth)

  ! now pre-compute the diffusion and reaction coefficients
  call Grid_getMinCellSize(dx)

  fl_R_over_s = fl_kpp_fact*4.0/fl_b/dx
  fl_kappa_over_s = fl_b*dx/16.0

  fl_width = fl_b*dx

  ! get the guardcell masking needs of the flamespeed computation
  fl_gcMask(:) = .false.
  fl_gcDoEos = .false.
  call fl_fsGcMask(fl_gcMask,fl_gcDoEos)
  ! and add flame progress variable
  fl_gcMask(FLAM_MSCALAR) = .true.
  ! and that is a mass scalar so we need density for interpolation at refinement boundaries
  fl_gcMask(DENS_VAR) = .true.

  fl_gcDoLogMask = .true.
  

  if(fl_useFlame) then 
  ! call effects init first in case flamespeed needs some of its info
     call fl_effInit()
     call fl_fsInit()
  end if

  return
end subroutine Flame_init
