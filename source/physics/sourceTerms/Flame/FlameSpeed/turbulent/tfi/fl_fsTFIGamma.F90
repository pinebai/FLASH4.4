!!****if* source/physics/sourceTerms/Flame/FlameSpeed/turbulent/tfi/fl_fsTFIGamma
!!
!! NAME
!!
!!  fl_fsTFIGamma
!!
!! SYNOPSIS
!!
!!  call fl_fsTFIGamma(real(out) :: g,
!!                     real(in) :: up_over_s,
!!                     real(in) :: de_over_dl)
!!
!! DESCRIPTION
!!
!! Aaron Jackson 2010
!!
!! This subroutine calculates gamma.
!!
!! ARGUMENTS
!!
!!   g : 
!!
!!   up_over_s : 
!!
!!   de_over_dl : 
!!
!!
!!
!!***


subroutine fl_fsTFIGamma(G, up_over_s, de_over_dl)

  implicit none

  real, intent(out) :: G
  real, intent(in) :: up_over_s, de_over_dl

  ! calculate Gamma, eq. 30
  G = 0.75e0 * exp( -1.2e0 / up_over_s**0.3 ) * de_over_dl**(2.0/3.0)

  return
end subroutine fl_fsTFIGamma
