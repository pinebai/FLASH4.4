!!****if* source/physics/sourceTerms/Flame/FlameSpeed/Constant/fl_flameSpeed
!!
!! NAME
!!
!!  fl_flameSpeed
!!
!! SYNOPSIS
!!
!!  call fl_flameSpeed(real, dimension(:,:,:,:),POINTER_INTENT_IN  :: solndata,
!!                     real, dimension(:,:,:)(out) :: flamespeed,
!!                     integer(in) :: blockid,
!!                     integer(in) :: nlayers)
!!
!! DESCRIPTION
!!
!! Dean Townsley 2008
!!
!! ARGUMENTS
!!
!!   solndata : solution data
!!
!!   flamespeed : flame speed
!!
!!   blockid : ID of block in current processor
!!
!!   nlayers : number of layers
!!
!!
!!
!!***


subroutine fl_flameSpeed( solnData, flamespeed, blockID, nlayers)

#include "Flash.h"
#include "constants.h"
#include "FortranLangFeatures.fh"

  use fl_fsData, only : fl_fsConstFlameSpeed, fl_fsConstFlameWidth
  implicit none
  real, dimension(:,:,:,:),POINTER_INTENT_IN :: solnData
  real, dimension(:,:,:),intent(out) :: flamespeed
  integer, intent(in) :: blockID, nlayers

  flamespeed(:,:,:) = fl_fsConstFlameSpeed
#ifdef FSPD_VAR
  solndata(FSPD_VAR,:,:,:) = fl_fsConstFlameSpeed
#endif

end subroutine
