!!****if* source/physics/Eos/EosMain/Gamma/eos_initGamma
!!
!! NAME
!!
!!  eos_initGamma
!!
!! SYNOPSIS
!!  
!!  call  eos_initGamma()
!!                 
!!
!! DESCRIPTION
!!
!!  This routine does ideal gamma law specific initialization
!!
!!  ARGUMENTS
!!
!!
!!***

#include "Eos.h"
subroutine eos_initGamma()

  use Eos_data, ONLY : eos_type
  use Eos_data, ONLY : eos_gamma
  use eos_idealGammaData, ONLY : eos_gammam1

  implicit none

  eos_type=EOS_GAM

  eos_gammam1 = 1.0/(eos_gamma-1.0)

  return
end subroutine eos_initGamma
