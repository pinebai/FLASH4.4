!!****if* source/physics/sourceTerms/Flame/FlameSpeed/turbulent/fl_fsTFIGcMask
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
!! Stub
!! Aaron Jackson 2010
!!
!! ARGUMENTS
!!
!!   fl_gcmask : flame guardcell mask
!!
!!   fl_gcdoeos : do eos check
!!
!!
!!
!!***

! See fl_fsTFIInterface for a description of subroutines
!
!

#include "Flash.h"
subroutine fl_fsTFIGcMask(fl_gcMask,fl_gcDoEos)

  implicit none

  logical, dimension(:), intent(inout) ::  fl_gcMask
  logical, intent(inout) :: fl_gcDoEos

  return

end subroutine
