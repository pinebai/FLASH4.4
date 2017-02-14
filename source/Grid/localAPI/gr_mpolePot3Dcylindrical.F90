!!****if* source/Grid/localAPI/gr_mpolePot3Dcylindrical
!!
!! NAME
!!
!!  gr_mpolePot3Dcylindrical
!!
!! SYNOPSIS
!!
!!  gr_mpolePot3Dcylindrical  (integer, intent(in) :: ipotvar)
!!
!! DESCRIPTION
!!
!!  Computes the potential field for a three-dimensional cylindrical geometry
!!  using the mass moments already calculated. On output the variable
!!  indexed by ipotvar contains the potential. The calculations are
!!  entirely local to each processor, since each processor has a local
!!  copy of the moments.
!!
!! ARGUMENTS
!!
!!  ipotvar : index to variable containing the potential
!!
!! SIDE EFFECTS
!!
!!  Insignificant side effects:
!!  (the initial and final values in these arrays do not matter.  The
!!  arrays exist at module scope so that they only need to be
!!  allocated once - they can be viewed as a persistent scratch space.
!!  The arrays are completely reinitialized in each subroutine which
!!  accesses the array values.)
!!   * gr_mpoleDampingC
!!   * gr_mpoleDampingI
!!   * gr_mpoleDampingR
!!   * gr_mpoleI
!!   * gr_mpoleR
!!
!!***

!!REORDER(4): solnData

subroutine gr_mpolePot3Dcylindrical (ipotvar)

  implicit none
  
  integer, intent (in) :: ipotvar

  return
end subroutine gr_mpolePot3Dcylindrical
