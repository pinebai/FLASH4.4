!!****if* source/Grid/localAPI/gr_mpolePot2Dcylindrical
!!
!! NAME
!!
!!  gr_mpolePot2Dcylindrical
!!
!! SYNOPSIS
!!
!!  gr_mpolePot2Dcylindrical  (integer, intent(in) :: ipotvar)
!!
!! DESCRIPTION
!!
!!  Computes the potential field for a two-dimensional cylindrical geometry
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

subroutine gr_mpolePot2Dcylindrical (ipotvar)

  implicit none
  
  integer, intent (in) :: ipotvar

  return
end subroutine gr_mpolePot2Dcylindrical

