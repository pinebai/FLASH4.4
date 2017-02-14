!!****if* source/Grid/localAPI/gr_mpoleCen1Dspherical
!!
!! NAME
!!
!!  gr_mpoleCM1Dspherical
!!
!! SYNOPSIS
!!
!!  gr_mpoleCen1Dspherical (integer, intent(in) :: idensvar)
!!
!! DESCRIPTION
!!
!!  Computes all data related to the center of multipole expansion for 1D spherical
!!  geometry. For a 1D spherical problem the center of expansion for the multipoles
!!  is always at the radial domain origin and does not need to be computed. The
!!  following is computed here:
!!
!!                  1) total mass (aborts, if <= 0)
!!                  2) the 'atomic' inner zone length (and its inverse)
!!
!! ARGUMENTS
!!
!!  idensvar : the index of the density variable
!!
!!***

!!REORDER(4): solnData

subroutine gr_mpoleCen1Dspherical (idensvar)

  implicit none
  
  integer, intent (in) :: idensvar

  return
end subroutine gr_mpoleCen1Dspherical
