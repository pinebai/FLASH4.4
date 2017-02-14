!!****if* source/physics/sourceTerms/Flame/FlameSpeed/LaminarOnly/fl_fsGcMask
!!
!! NAME
!!
!!  fl_fsGcMask
!!
!! SYNOPSIS
!!
!!  call fl_fsGcMask(logical, dimension(:)(inout) :: fl_gcmask,
!!                   logical(inout) :: fl_gcdoeos)
!!
!! DESCRIPTION
!!
!! Dean Townsley 2008
!!
!!
!! ARGUMENTS
!!
!!   fl_gcmask : GC mask
!!
!!   fl_gcdoeos : Check for Eos
!!
!!
!!
!!***

#include "Flash.h"
subroutine fl_fsGcMask(fl_gcMask,fl_gcDoEos)

  use fl_fsTFIInterface, ONLY : fl_fsTFIGcMask

  implicit none

  logical, dimension(:), intent(inout) ::  fl_gcMask
  logical, intent(inout) :: fl_gcDoEos

  ! pressure is needed to estimate unburned denisty
  ! which is in turn used to calculate laminar flame speed
  fl_gcMask(PRES_VAR) = .true.

  ! for now we'll set the gc filling code do the eos
  ! might consider doing it ourselves later to save work
  fl_gcDoEos = .true.

  call fl_fsTFIGcMask(fl_gcMask, fl_gcDoEos)

  return

end subroutine
