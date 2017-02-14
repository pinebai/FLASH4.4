!!****if* source/physics/sourceTerms/Flame/FlameSpeed/BuoyancyCompensation/fl_fsGcMask
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
!! Guardcell mask
!!
!! ARGUMENTS
!!
!!   fl_gcmask : Guardcell mask
!!
!!   fl_gcdoeos : Eos check
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

  ! gravitational potential is needed for Buoyancy compensation
  fl_gcMask(GPOT_VAR) = .true.
  ! pressure is needed to estimate unburned denisty
  ! which is in turn used to calculate laminar flame speed
  fl_gcMask(PRES_VAR) = .true.

  ! need initial abundance for flame speed and width and atwood number
#ifdef CI_MSCALAR
  fl_gcMask(CI_MSCALAR) = .true.
#endif
#ifdef NEI_MSCALAR
  fl_gcMask(NEI_MSCALAR) = .true.
#endif
  ! those are mass scalars so we need density for interppolation too
  fl_gcMask(DENS_VAR) = .true.

  ! the pressure in the guard cells is used to estimate the
  ! unburned density, so (1) it doesn't need to be
  ! uber-accurate since it is an estimate anyway and
  ! (2) the neighboring block is likely to be at the
  ! same refinement level near the flame anyway, in
  ! which case the eos call is not even necessary
  !  this defaults to false so we don't actually set it, this
  !  note is just here so we know why we don't
  !fl_gcDoEos = .false.

  call fl_fsTFIGcMask(fl_gcMask, fl_gcDoEos)

  return

end subroutine
