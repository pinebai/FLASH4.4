!!****if* source/physics/sourceTerms/Flame/FlameMain/Flame_getWidth
!!
!! NAME
!!
!!  Flame_getWidth
!!
!! SYNOPSIS
!!
!!  call Flame_getWidth(real, intent(OUT)  :: laminarwidth)
!!
!! DESCRIPTION
!!
!! Dean Townsley 2008
!!  approximate total width of the flame front
!!  more than about twice this far away progress
!!  variable can be initialized  to 0 or 1 for
!!  unburned and burned respectively
!!
!! ARGUMENTS
!!
!!   laminarwidth : 
!!
!!
!!
!!***

subroutine Flame_getWidth(laminarWidth)

  use Flame_data, ONLY : fl_width

  implicit none
  real, intent(OUT) :: laminarWidth

  laminarWidth=fl_width
  return

end subroutine Flame_getWidth
