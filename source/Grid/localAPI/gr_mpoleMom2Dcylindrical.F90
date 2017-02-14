!!****if* source/Grid/localAPI/gr_mpoleMom2Dcylindrical
!!
!! NAME
!!
!!  gr_mpoleMom2Dcylindrical
!!
!! SYNOPSIS
!!
!!  gr_mpoleMom2Dcylindrical (integer, intent(in) :: idensvar)
!!
!! DESCRIPTION
!!
!!  Compute the multipole moments of the density distribution in a
!!  two-dimensional cylindrical geometry, assuming the center of mass and
!!  the total mass have first been computed. On output, the Moment_R ()
!!  and Moment_I () arrays contain the mass moments over the regular
!!  and irregular solid harmonics.
!!
!! ARGUMENTS
!!
!!  idensvar -- the index of the density variable
!!
!!***

!!REORDER(4): solnData

subroutine gr_mpoleMom2Dcylindrical (idensvar)

  implicit none
  
  integer, intent(in) :: idensvar

  return
end subroutine gr_mpoleMom2Dcylindrical
