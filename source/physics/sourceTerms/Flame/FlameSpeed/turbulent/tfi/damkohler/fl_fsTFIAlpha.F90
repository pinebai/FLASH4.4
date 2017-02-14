!!****if* source/physics/sourceTerms/Flame/FlameSpeed/turbulent/tfi/damkohler/fl_fsTFIAlpha
!!
!! NAME
!!
!!  fl_fsTFIAlpha
!!
!! SYNOPSIS
!!
!!  call fl_fsTFIAlpha(real(out) :: alpha,
!!                     real(in) :: up,
!!                     real(in) :: de,
!!                     real(in) :: dl0)
!!
!! DESCRIPTION
!!
!! Aaron Jackson 2010
!!
!! This subroutine calculates the model constant alpha in the
!! Colin et al. (2000) TFI perscription
!!
!! ARGUMENTS
!!
!!   alpha : 
!!
!!   up : 
!!
!!   de : 
!!
!!   dl0 : 
!!
!!
!!
!!***


subroutine fl_fsTFIAlpha(alpha, up, de, dl0)

#include "Flash.h"

  use fl_fsTFIData, ONLY : fl_fsTFIBeta

  implicit none

  real, intent(out) :: alpha
  real, intent(in) :: up, de, dl0

  ! fl_fsBeta already includes 2 ln2 / 3 c_ms pre-factor
  ! this assumes we recover Damkohler scaling
  alpha = fl_fsTFIBeta * ( dl0 / de )**(2.0/3.0) 

  return

end subroutine fl_fsTFIAlpha
