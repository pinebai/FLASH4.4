!!****if* source/Grid/GridMain/paramesh/gr_xyzToBlockLevel
!!
!!  NAME
!!     gr_xyzToBlockLevel
!!
!!  SYNOPSIS
!!     call gr_xyzToBlockLevel(integer(in)  :: lev,
!!                              real  (in)  :: xyz(NDIM),
!!                             integer(OUT) :: ijk(NDIM))
!!
!!  DESCRIPTION
!!
!!    Maps physical domain coordinates to integer indices at a specific
!!    refinement level.  The valid range for an index is level dependent.
!!    Specifically, an index locating a block in the domain will be in the
!!    range [ 0, 2**(lev-1) ) if all relevant values of gr_nblockX,
!!    gr_nblockY,gr_nblockZ are 1; the range is increased appropriately
!!    for larger values of gr_nblockX,gr_nblockY,gr_nblockZ.
!!    Physical coordinates outside of the domain
!!    will map to indices outside the above range.
!!
!!  ARGUMENTS
!!    lev: (in) 1-based level we want i,j,k block coordinates on
!!    xyz: (in) point in the domain we want mapped to a block
!!    ijk: (out) 0-based integer coordinate of block for the level requested
!!
!!***
subroutine gr_xyzToBlockLevel(lev, xyz, ijk)
#include "Flash.h"
  use Grid_data, only: &
    gr_imin, gr_jmin, gr_kmin, &
    gr_imax, gr_jmax, gr_kmax, &
    gr_nblockX, gr_nblockY, gr_nblockZ
  implicit none
  
  integer, intent(in) :: lev
  real, intent(in) :: xyz(NDIM)
  integer, intent(out) :: ijk(NDIM)
  
  integer :: top(NDIM)
  real :: gmin(NDIM), gmax(NDIM)
  top = int(reshape((/gr_nblockX, gr_nblockY, gr_nblockZ/), (/NDIM/)))
  gmin = reshape((/gr_imin, gr_jmin, gr_kmin/), (/NDIM/))
  gmax = reshape((/gr_imax, gr_jmax, gr_kmax/), (/NDIM/))
  ijk = floor((xyz-gmin)/(gmax-gmin)*ishft(top,lev-1))
end subroutine gr_xyzToBlockLevel
