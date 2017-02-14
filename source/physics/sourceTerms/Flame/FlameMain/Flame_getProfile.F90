!!****if* source/physics/sourceTerms/Flame/FlameMain/Flame_getProfile
!!
!! NAME
!!
!!  Flame_getProfile
!!
!! SYNOPSIS
!!
!!  call Flame_getProfile(real(in) :: x,
!!                        real(out) :: f)
!!
!! DESCRIPTION
!!
!!   get value of progress variable a distance x from center of
!!   flame front in the steady state propagating flame (used toDean Townsley 2008
!!   initialize data on mesh).
!!
!! ARGUMENTS
!!   x : x is defined such that positive x is in the direction of propagation
!!
!!   f :  f = 0.5 at x = 0
!!
!!
!!
!!***

subroutine Flame_getProfile(x, f)

  use Flame_data, ONLY : fl_width, fl_initProfileAdjustWidth

  implicit none
  real, intent(in)  :: x
  real, intent(out) :: f

  ! This is an approximate profile form based on widths determined in
  ! Vladimirova et al.  Here the "width" is approximately the
  ! distance between where phi=0.12 and 0.88.  Over twice width, phi
  ! goes from 0.02 to 0.98.
  ! The fl_initProfileAdjustmentWidth is to allow compatibility with
  ! slight variations if necessary

  ! tanh is slow, but no need to be super efficient here since this
  ! should only be called at init
  f = 0.5 * (1.0 - tanh(x/fl_width/0.5/fl_initProfileAdjustWidth))

  return

end subroutine Flame_getProfile
