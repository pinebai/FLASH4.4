!!****if* source/physics/RadTrans/RadTransMain/MGD/RadTrans_mgdSetBound
!!
!!  NAME 
!!
!!  RadTrans_mgdSetBound
!!
!!  SYNOPSIS
!!
!!  call RadTrans_mgdSetBound( integer(IN) :: g,
!!                             real(IN) :: b )
!!
!!  DESCRIPTION 
!!      This subroutine is used to set a particular radiation
!!      energy group boundary for MGD.
!!      DEV : What does this mean?
!!
!! ARGUMENTS
!!
!!      g : The boundary number, group g is bounded by g and g+1
!!      b : The boundary energy [ergs]
!! 
!!***
subroutine RadTrans_mgdSetBound(g, b)
  use rt_data, ONLY: rt_mgdBounds
  implicit none

  integer, intent(in) :: g
  real,    intent(in) :: b
  rt_mgdBounds(g) = b
  return

end subroutine RadTrans_mgdSetBound
