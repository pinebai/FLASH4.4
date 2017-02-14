!!****f* source/physics/RadTrans/RadTrans_mgdSetBound
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
!!
!! ARGUMENTS
!!
!!      g : The boundary number, group g is bounded by g and g+1
!!      b : The boundary energy [ergs]
!! 
!!***
subroutine RadTrans_mgdSetBound(g, b)
  implicit none

  integer, intent(in) :: g
  real,    intent(in) :: b
  ! Stub implementation
  return

end subroutine RadTrans_mgdSetBound
