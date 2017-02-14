!!****if* source/physics/sourceTerms/Flame/FlameSpeed/turbulent/tfi/fl_fsTFIEnhance
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


subroutine fl_fsTFIEnhance(E, up, s, de, dl0, de_over_dl1, E_lim)

  use fl_fsTFIInterface, ONLY : fl_fsTFIAlpha, fl_fsTFIGamma

  implicit none

  real, intent(out) :: E
  real, intent(in) :: up, s, de, dl0, de_over_dl1
  real, intent(out), optional :: E_lim

  real :: up_over_s, de_over_dl0, a, G0, G1, W0, W1

  ! use negative to represent no quenching limit
  if (present(E_lim)) E_lim = -1.0e0

  up_over_s = up / s
  de_over_dl0 = de / dl0

  ! get alpha
  call fl_fsTFIAlpha(a, up, de, dl0)

  ! calculate Gamma, eq. 30
  call fl_fsTFIGamma(G0, up_over_s, de_over_dl0)
  call fl_fsTFIGamma(G1, up_over_s, de_over_dl1)

  ! calculate wrinkling, eq. 32
  W0 = 1.0e0 + a * G0 * up_over_s
  W1 = 1.0e0 + a * G1 * up_over_s

  ! calculate enhancement
  E = W0 / W1

  return
end subroutine fl_fsTFIEnhance
