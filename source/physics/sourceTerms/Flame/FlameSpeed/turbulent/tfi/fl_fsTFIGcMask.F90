!!****if* source/physics/sourceTerms/Flame/FlameSpeed/turbulent/tfi/fl_fsTFIGcMask
!!
!! NAME
!!
!!  fl_fsTFIGcMask
!!
!! SYNOPSIS
!!
!!  call fl_fsTFIGcMask(logical, dimension(:)(inout) :: fl_gcmask,
!!                      logical(inout) :: fl_gcdoeos)
!!
!! DESCRIPTION
!!
!! Aaron Jackson 2010
!!
!!
!! ARGUMENTS
!!
!!   fl_gcmask : guardcell mask
!!
!!   fl_gcdoeos : check for eos
!!
!!
!!
!!***


#include "Flash.h"
subroutine fl_fsTFIGcMask(fl_gcMask,fl_gcDoEos)

  implicit none

  logical, dimension(:), intent(inout) ::  fl_gcMask
  logical, intent(inout) :: fl_gcDoEos

  ! turbulent strength is needed to calculate the TFI-enhanced flame speed
  fl_gcMask(TURB_VAR) = .true.

  return

end subroutine
