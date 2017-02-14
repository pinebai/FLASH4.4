!!****f* source/physics/materialProperties/MagneticResistivity/MagneticResistivity_fullState
!!
!! NAME
!!  MagneticResistivity_fullState
!!
!! SYNOPSIS
!!  call MagneticResistivity_fullState(real(in)    :: solnVec(NUNK_VARS),
!!                                     real(out)   :: resPar,
!!                            OPTIONAL,real(out)   :: resPerp)
!!
!! DESCRIPTION
!!
!! Computes the Spitzer electron Magnetic Resistivity for all materials,
!! including those with Z > 1. The specific equations used here all
!! come from "The Physics of Inertial Fusion" by Atzeni.
!!
!!  The stub implementation returns resPar = 0
!!
!! ARGUMENTS
!!
!!   solnVec  :   solution state, a vector from UNK with all variables
!!   resPar   :   parallel component of Magnetic Resistivity
!!   resPerp :    perpendicular component of Magnetic Resistivity
!!
!!***

#include "Flash.h"

subroutine MagneticResistivity_fullState(solnVec,resPar, resPerp)
  implicit none

  real, intent(in)  :: solnVec(NUNK_VARS)
  real, intent(out) :: resPar
  real, OPTIONAL, intent(out) :: resPerp

  ! Stub
  resPar  = 0.0
  if (present(resPerp)) resPerp = 0.0
end subroutine MagneticResistivity_fullState
