!!****f* source/physics/RadTrans/RadTrans_mgdGetBound
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
  implicit none

  integer, intent(in) :: g
  real,    intent(out) :: b

  ! Stub implementation
  b = 0.0
  return

end subroutine RadTrans_mgdGetBound
