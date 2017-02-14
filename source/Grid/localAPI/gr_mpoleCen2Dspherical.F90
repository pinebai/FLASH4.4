!!****if* source/Grid/localAPI/gr_mpoleCen2Dspherical
!!
!! NAME
!!
!!  gr_mpoleCen2Dspherical
!!
!! SYNOPSIS
!!
!!  gr_mpoleCen2Dspherical (integer, intent(in) :: idensvar)
!!
!! DESCRIPTION
!!
!!  Computes all data related to the center of multipole expansion for 2D spherical
!!  geometry. It computes the center of expansion for the multipoles for 2D spherical
!!  geometries. The center is calculated using the position weighted by the square
!!  density:
!!
!!
!!                          integral (r * rho * rho  dr)
!!                Cen (z) = ----------------------------
!!                            integral (rho * rho  dr)
!!
!!
!!  which, due to uniform density in each cell, becomes:
!!
!!
!!                   sum cells (cell center r * cell mass * cell rho)
!!         Cen (z) = ------------------------------------------------
!!                           sum cells (cell mass * cell rho)
!!
!!
!!  After the initial Cen (z) has been determined, it is placed on the nearest
!!  nearest cell corner. The following is computed here:
!!
!!                  1) multipole expansion center (placed on nearest cell corner)
!!                  2) total mass (aborts, if <= 0)
!!                  3) the 'atomic' inner zone length (and its inverse)
!!
!! ARGUMENTS
!!
!!  idensvar : the index of the density variable
!!
!! NOTES
!!
!!  gr_mpoleZcenter denotes the location of the center of multipole expansion
!!  in the 2D cartesian version of the 2D spherical framework (R,theta -> x,z).
!!  Due to symmetry, there is no need to evaluate and store a x-coordinate
!!  component, which is always equal to zero.
!!
!!***

!!REORDER(4): solnData

subroutine gr_mpoleCen2Dspherical (idensvar)

  implicit none
  
  integer, intent (in) :: idensvar

  return
end subroutine gr_mpoleCen2Dspherical
