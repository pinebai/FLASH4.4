!!****if* source/physics/sourceTerms/Flame/FlameSpeed/turbulent/tfi/pocheau/fl_fsTFIAlpha
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

  implicit none

  real, intent(out) :: alpha
  real, intent(in) :: up, de, dl0

  return

end subroutine fl_fsTFIAlpha
