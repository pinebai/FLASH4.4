!!****if* source/physics/sourceTerms/Flame/FlameSpeed/Constant/TFI/fl_fsGcMask
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


subroutine fl_fsGcMask(fl_gcMask,fl_gcDoEos)

  use fl_fsTFIInterface, ONLY : fl_fsTFIGcMask

  implicit none

  logical, dimension(:), intent(inout) ::  fl_gcMask
  logical, intent(inout) :: fl_gcDoEos

  ! this flame speed doesn't need anything in the guardcells
  ! so we leave both alone, as they default to false

  call fl_fsTFIGcMask(fl_gcMask, fl_gcDoEos)

  return

end subroutine
