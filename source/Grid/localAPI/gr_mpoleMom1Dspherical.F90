!!****if* source/Grid/localAPI/gr_mpoleMom1Dspherical
!!
!! NAME
!!
!!  gr_mpoleMom1Dspherical
!!
!! SYNOPSIS
!!
!!  gr_mpoleMom1Dspherical (integer, intent(in) :: idensvar)
!!
!! DESCRIPTION
!!
!!  Compute the multipole moments of the density distribution for
!!  a one-dimensional spherical geometry, assuming the center of mass
!!  and the total mass have first been computed.  On output, the
!!  Moment_R () and Moment_I () arrays contain the mass moments
!!  over the regular and irregular solid harmonics.
!!
!! ARGUMENTS
!!
!!  idensvar -- the index of the density variable
!!
!!***

!!REORDER(4): solnData

subroutine gr_mpoleMom1Dspherical (idensvar)

  implicit none
  
  integer, intent (in) :: idensvar

  return
end subroutine gr_mpoleMom1Dspherical
