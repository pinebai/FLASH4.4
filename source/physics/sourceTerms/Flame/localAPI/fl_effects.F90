!!****if* source/physics/sourceTerms/Flame/localAPI/fl_effects
!!
!! NAME
!!
!!  fl_effects
!!
!! SYNOPSIS
!!
!!  call fl_effects(real, dimension(:,:,:,:),POINTER_INTENT_IN  :: solndata,
!!                  real,dimension(:,:,:)(in) :: flamdot,
!!                  real(in) :: dt,
!!                  integer(in) :: blockid)
!!
!! DESCRIPTION
!! Dean Townsley 2008
!!
!! this is a stub for flame effects that don't do anything locally
!!
!! ARGUMENTS
!!
!!   solndata : solution data 
!!
!!   flamdot : 
!!
!!   dt : timestep
!!
!!   blockid : ID of block in current processor
!!
!!
!!
!!***


#include "FortranLangFeatures.fh"

subroutine fl_effects( solnData, flamdot, dt, blockID)

  implicit none

  real, dimension(:,:,:,:),POINTER_INTENT_IN  :: solnData
  real,dimension(:,:,:), intent(in)     :: flamdot
  real,intent(in)                       :: dt
  integer, intent(in)                   :: blockID

  return
end subroutine
