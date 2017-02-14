!!****f* source/physics/sourceTerms/Flame/Flame_getProfile
!!
!! NAME
!!
!!  Flame_getProfile
!!
!! SYNOPSIS
!!
!!  call Flame_getProfile (real,intent(in) :: x, real,intent(out) :: f)
!!
!! DESCRIPTION
!!
!!  Get value of progress variable a distance x from center of
!!  flame front in the steady state propagating flame (used to
!!  initialize data on mesh).
!!  x is defined such that positive x is in the direction of
!!  propagation and f = 0.5 at x = 0
!!
!! ARGUMENTS
!!
!!   x - distance from flame center.
!!   f - flame progress variable value returned.
!!
!! SEE ALSO
!!
!!  see Flame_interface.F90 for possible updates
!!
!!***

! this is a stub for when the Flame Unit is not included
!
! Dean Townsley 2008
!
subroutine Flame_getProfile(x, f)
  implicit none
  real, intent(in)  :: x
  real, intent(out) :: f

  f=0.0
  return

end subroutine Flame_getProfile
