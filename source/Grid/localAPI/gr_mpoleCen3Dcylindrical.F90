!!****if* source/Grid/localAPI/gr_mpoleCen3Dcylindrical
!!
!! NAME
!!
!!  gr_mpoleCen3Dcylindrical
!!
!! SYNOPSIS
!!
!!  gr_mpoleCen3Dcylindrical (integer, intent (in) :: idensvar)
!!
!! DESCRIPTION
!!
!!  Computes all data related to the center of multipole expansion for 3D cylindrical
!!  geometry. It computes the center of expansion for the multipoles for 3D cylindrical
!!  geometries. The center is calculated using the position weighted by the square
!!  density:
!!
!!
!!                            integral (r * rho * rho  dr)
!!              Cen (x,y,z) = ----------------------------
!!                              integral (rho * rho  dr)
!!
!!
!!  which, due to uniform density in each cell, becomes:
!!
!!
!!                     sum cells (cell center r * cell mass * cell rho)
!!       Cen (x,y,z) = ------------------------------------------------
!!                             sum cells (cell mass * cell rho)
!!
!!
!!  After the initial Cen (x,y,z) has been determined, it is placed on the
!!  the nearest cell corner. The following is computed here:
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
!!  gr_mpoleXcenter, gr_mpoleYcenter and gr_mpoleZcenter denote the location of
!!  the center of multipole expansion in the 3D cartesian framework (going from
!!  R,phi,z -> x,y,z).
!!
!!***

!!REORDER(4): solnData

subroutine gr_mpoleCen3Dcylindrical (idensvar)

  implicit none
  
  integer, intent (in) :: idensvar

  return
end subroutine gr_mpoleCen3Dcylindrical
