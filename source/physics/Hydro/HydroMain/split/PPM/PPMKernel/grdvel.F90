!!****if* source/physics/Hydro/HydroMain/split/PPM/PPMKernel/grdvel
!!
!! NAME
!!
!!  grdvel
!!
!! SYNOPSIS
!!
!!  grdvel(integer(IN):: blockID)
!!
!! DESCRIPTION
!!
!!  Compute grid velocity for moving-grid artificial viscosity
!!  method in PPM.
!!
!!  Currently this routine appears to assume a sweep in the
!!  x-direction for a 2D grid.  1D, 3D, and 2D grids with sweeps
!!  in the y-direction don't appear to be supported.
!!
!!  I really don't see how this function could possible work -- the grid velocity
!!  is going to be block dependent, so it is going to be different at different
!!  points in the computational domain.  For now, this should not be used, which
!!  is accomplished by setting vgrid to 0 in flash.par.  Eventually, this should
!!  be ripped out of the code.
!!
!! 
!! ARGUMENTS
!! 
!! blockID : Local block ID.
!!
!! PARAMETERS
!!
!!
!!***

subroutine grdvel (blockID)

#include "Flash.h"

  use Hydro_data, ONLY : hy_vgrid

  implicit none

  integer, intent(in) :: blockID

  real :: vmax

  !integer, save :: ivelx, ively, iugrid
  !integer, save :: iXcoord, iYcoord, iZcoord

  integer :: i, j, k

!!$  real, DIMENSION(IHI_GC) :: ugridx
!!$  real, DIMENSION(JHI_GC) :: ugridy
!!$  real, DIMENSION(KHI_GC) :: ugridz
!!$
!!$! we will access the solution data through a pointer return by the 
!!$! accessor method
!!$  real, DIMENSION(:,:,:,:), POINTER :: solnData
!!$  logical :: gcell = .true.
!!$
!!$
!!$!  iugrid = dBaseKeyNumber('ugrid')
!!$     
!!$
!!$  !call get_parm_from_context(GLOBAL_PARM_CONTEXT, 'vgrid', vgrid)
!!$
!!$
!!$  !!solnData => dBaseGetDataPtrSingleBlock(blockID, GC)
!!$  call Grid_getBlkPtr(blockID, solnData)
!!$
!!$!===============================================================================
!!$
!!$! start by computing the maximum grid velocity -- this is block dependent ????
!!$
!!$! this seems to assume that the entire computation domain is represented by a
!!$! single patch -- this goes back to the original serial PROMETHEUS code.  The
!!$! following loop was supposed to look at the velocity at the edge of the 
!!$! computation domain and find the maximum, to determine how much to expand the
!!$! grid by.  Now it looks at each block, even if they are in the middle of the 
!!$! domain.  Furthermore, the vmax is not shared across processors or blocks, so
!!$! every block will have a different grid velocity -- this cannot be good.
!!$
!!$  ugridx(:) = 0.e0
!!$  ugridy(:) = 0.e0
!!$  ugridz(:) = 0.e0
!!$
!!$  vmax = 0.e0
!!$
!!$  do j = JLO, JHI
!!$     vmax = max (vmax, solnData(VELX_VAR,NGUARD+NXB,j,1))
!!$  end do
!!$
!!$  vmax = vmax * hy_vgrid
!!$
!!$! x grid velocity
!!$  do i = ILO, IHI
!!$     ugridx(i) = vmax * (i - 0.5e00) / (NGUARD+NXB - 0.5e00)
!!$  end do
!!$
!!$
!!$  !!!call dBasePutCoords(iugrid, iXcoord, blockID, ugridx)
!!$
!!$
!!$! y grid velocity
!!$  do j = JLO, JHI
!!$     ugridy(j) = 0.e0
!!$  end do
!!$
!!$
!!$  !!call dBasePutCoords(iugrid, iYcoord, blockID, ugridy)
!!$
!!$
!!$! z grid velocity
!!$  do k = KLO, KHI
!!$     ugridz(k) = 0.e0
!!$  end do
!!$
!!$  !!!call dBasePutCoords(iugrid, iZcoord, blockID, ugridz)
!!$
!!$
!!$
!!$  !!call dBaseReleaseDataPtrSingleBlock(blockID, solnData)
!!$  call Grid_releaseBlkPtr(blockID,solnData)

  return
end subroutine grdvel
