!!****if* source/physics/RadTrans/RadTransMain/MGD/RadTrans_mgdGetBound
!!
!!  NAME 
!!
!!  RadTrans_mgdGetBound
!!
!!  SYNOPSIS
!!
!!  call RadTrans_mgdGetBound( integer(IN) :: g,
!!                             real(OUT) :: b )
!!
!!  DESCRIPTION 
!!      This subroutine is used to access a particular radiation
!!      energy group boundary for MGD.
!!
!! ARGUMENTS
!!
!!      g : The boundary number, group g is bounded by g and g+1
!!      b : The boundary energy [ergs]
!! 
!!***
subroutine RadTrans_mgdGetBound(g, b)
  use rt_data, ONLY: rt_mgdBounds
  implicit none

  integer, intent(in) :: g
  real,    intent(out) :: b
  b = rt_mgdBounds(g)
  return

end subroutine RadTrans_mgdGetBound
