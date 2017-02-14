!!****if* source/physics/sourceTerms/Turb/TurbMain/Turb_init
!!
!! NAME
!!
!!  Turb_init
!!
!! SYNOPSIS
!!
!!  call Turb_init()
!!
!! DESCRIPTION
!!
!! Dean Townsley 2008
!!
!! ARGUMENTS
!!
!!   No arguments
!!
!!
!!
!!***


#include "Flash.h"
subroutine Turb_init()

  use Driver_interface, only : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Turb_data

  implicit none

  ! get runtime parameters
  call RuntimeParameters_get("useTurb", turb_useTurb)
  call RuntimeParameters_get("turb_stepSize", turb_stepSize)
  call RuntimeParameters_get("turb_c2", turb_c2)

  if (2*turb_stepSize > NGUARD) &
     call Driver_abortFlash("stencil size of turb operator > number of guard cells")

  ! check if both operators can be performed without guard cell filling
  turb_fillGC = (4*turb_stepSize > NGUARD)

  return
end subroutine Turb_init
