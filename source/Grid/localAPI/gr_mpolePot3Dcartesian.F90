!!****if* source/Grid/localAPI/gr_mpolePot3Dcartesian
!!
!! NAME
!!
!!  gr_mpolePot3Dcartesian
!!
!! SYNOPSIS
!!
!!  gr_mpolePot3Dcartesian  (integer, intent(in) :: ipotvar)
!!
!! DESCRIPTION
!!
!!  Computes the potential field for a three-dimensional cartesian geometry
!!  using the mass moments already calculated. On output the variable
!!  indexed by ipotvar contains the potential. The calculations are
!!  entirely local to each processor, since each processor has a local
!!  copy of the moments.
!!
!! ARGUMENTS
!!
!!  ipotvar : index to variable containing the potential
!!
!!***

!!REORDER(4): solnData

subroutine gr_mpolePot3Dcartesian (ipotvar)

  implicit none
  
  integer, intent (in) :: ipotvar

  return
end subroutine gr_mpolePot3Dcartesian
