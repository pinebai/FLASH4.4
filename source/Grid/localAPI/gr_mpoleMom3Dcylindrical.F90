!!****if* source/Grid/localAPI/gr_mpoleMom3Dcylindrical
!!
!! NAME
!!
!!  gr_mpoleMom3Dcylindrical
!!
!! SYNOPSIS
!!
!!  gr_mpoleMom3Dcylindrical (integer, intent(in) :: idensvar)
!!
!! DESCRIPTION
!!
!!  Compute the multipole moments of the density distribution in three-dimensional
!!  cylindrical geometry, assuming the center of mass and the total mass have first
!!  been computed. On output, the gr_mpoleMomentR and gr_mpoleMomentI arrays contain
!!  the mass moments over the regular and irregular solid harmonics.
!!
!! ARGUMENTS
!!
!!  idensvar : the index of the density variable
!!
!! SIDE EFFECTS
!!
!!  Significant side effects:
!!   * gr_mpoleMomentI
!!   * gr_mpoleMomentR
!!
!!  Insignificant side effects:
!!  (the initial and final values in these arrays do not matter.  The
!!  arrays exist at module scope so that they only need to be
!!  allocated once - they can be viewed as a persistent scratch space.
!!  The arrays are completely reinitialized in each subroutine which
!!  accesses the array values.)
!!   * gr_mpoleI
!!   * gr_mpoleR
!!
!!***

!!REORDER(4): solnData

subroutine gr_mpoleMom3Dcylindrical (idensvar)

  implicit none
  
  integer, intent (in) :: idensvar

  return
end subroutine gr_mpoleMom3Dcylindrical
