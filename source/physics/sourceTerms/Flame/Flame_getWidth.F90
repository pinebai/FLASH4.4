!!****f* source/physics/sourceTerms/Flame/Flame_getWidth
!!
!! NAME
!!
!!  Flame_getWidth
!!
!! SYNOPSIS
!!
!!  call Flame_getWidth ( real, intent(out) :: laminarWidth )
!!
!! DESCRIPTION
!!
!!  Gets the laminar width of artificial reaction-diffusion front.
!!
!!  Approximate total width of the flame front.
!!  More than about twice this far away progress
!!  variable con be initialized  to 0 or 1 for
!!  unburned and burned respectively
!!
!! ARGUMENTS
!!
!!  laminarWidth  -  the returned laminar width.
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
subroutine Flame_getWidth(laminarWidth)

  implicit none
  real, intent(OUT) :: laminarWidth

  laminarWidth=0.0
  return

end subroutine Flame_getWidth
