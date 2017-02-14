!!****if* source/physics/sourceTerms/Flame/FlameSpeed/turbulent/tfi/charlette/fl_fsTFIEnhance
!!
!! NAME
!!
!!  fl_fsTFIEnhance
!!
!! SYNOPSIS
!!
!!  call fl_fsTFIEnhance(real(out) :: e,
!!                       real(in) :: up,
!!                       real(in) :: s,
!!                       real(in) :: de,
!!                       real(in) :: dl0,
!!                       real(in) :: de_over_dl1,
!!                       real(out), optional  :: e_lim)
!!
!! DESCRIPTION
!!
!! Aaron Jackson 2010
!!
!! This subroutine calculates the enhancement factor to the flame speed.
!!
!! ARGUMENTS
!!
!!   e : 
!!
!!   up : 
!!
!!   s : 
!!
!!   de : 
!!
!!   dl0 : 
!!
!!   de_over_dl1 : 
!!
!!   e_lim : 
!!
!!
!!
!!***


#include "constants.h"

subroutine fl_fsTFIEnhance(E, up, s, de, dl0, de_over_dl1, E_lim)

  use fl_fsTFIInterface, ONLY : fl_fsTFIGamma
  use fl_fsTFIData, ONLY : fl_fsTFIBeta

  implicit none

  real, intent(out) :: E
  real, intent(in) :: up, s, de, dl0, de_over_dl1
  real, intent(out), optional :: E_lim
  real, parameter :: Ck = 1.5

  real :: up_over_s, de_over_dl0, fu, fD, fRe, a, b, G

  up_over_s = up / s
  de_over_dl0 = de / dl0

  call fl_fsTFIGamma(G, up_over_s, de_over_dl0)

  ! calculate enhancement
  if (present(E_lim)) E_lim = (1.0 + de_over_dl0)**fl_fsTFIBeta
  E = (1.0 + min(de_over_dl0,G*up_over_s))**fl_fsTFIBeta

  return
end subroutine fl_fsTFIEnhance
