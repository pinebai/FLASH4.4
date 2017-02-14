!!****if* source/Grid/GridSolvers/MultigridMC/gr_mgNorm
!!
!! NAME
!!
!!  gr_mgNorm
!!
!! SYNOPSIS
!!
!!  call gr_mgNorm(integer(in) :: level,
!!                 integer(in) :: ivar,
!!                 real(inout) :: norm,
!!                 integer(in) :: leaf_only)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   level : 
!!
!!   ivar : 
!!
!!   norm : 
!!
!!   leaf_only : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

!*******************************************************************************

!  Routine:     mg_norm()

!  Description: Compute the L2 norm of a multigrid variable on a particular
!               level.  Assume data on this level is good (e.g., because of
!               a call to mg_restrict()).  If level == 0, compute norm on
!               all levels.  If leaf_only /= 0, only compute norm on leaf
!               nodes.  If ivar == -1, compute norm on the work array.


subroutine gr_mgNorm (level, ivar, norm, leaf_only)

!===============================================================================
#include "Flash.h"

use gr_mgData, ONLY: nodetype_save, ili, jli, kli, iui, jui, kui

  use Grid_interface,    ONLY : Grid_getBlkPtr,       &
                                Grid_releaseBlkPtr,   &
                                Grid_getBlkPhysicalSize

use workspace, ONLY: work
use tree, only : nodetype,newchild,lrefine

implicit none
#include "Multigrid.h"
#include "constants.h"
include "mpif.h"

integer, intent(in) :: level, ivar, leaf_only
real, intent(inout) :: norm

integer :: lb, i, j, k, ierr, lnblocks
real    :: lvol, vol, lsum, bsum, sum
real    :: nbinv, cvol, bvol
logical :: include_in_sum

integer, parameter :: MAXDIMS = 3
real, dimension(MAXDIMS) :: size

real, pointer, dimension(:,:,:,:), save :: unk

integer, parameter :: nxb = NXB
integer, parameter :: nyb = NYB
integer, parameter :: nzb = NZB
integer, parameter :: nguard = NGUARD
integer, parameter :: ndim = NDIM
integer, parameter :: k2d = K2D
integer, parameter :: k3d = K3D

!===============================================================================

lvol = 0.
lsum = 0.
nbinv = 1. / (real(nxb)*real(nyb)*real(nzb))

call Grid_getLocalNumBlks(lnblocks)

do lb = 1, lnblocks
  include_in_sum = (lrefine(lb) == level) .or. (level == 0)
  if (leaf_only /= 0) & 
    include_in_sum = include_in_sum .and. (nodetype_save(lb) == 1)

  if (include_in_sum) then

    ! Get BlockSize:
    call Grid_getBlkPhysicalSize(lb,size)

    bvol = size(1)
    if (ndim >= 2) bvol = bvol * size(2)
    if (ndim == 3) bvol = bvol * size(3)
    cvol = bvol * nbinv
    lvol = lvol + bvol
    bsum = 0.
    if (ivar >= 0) then
    ! Point to blocks center vars:
    call Grid_getBlkPtr(lb,unk,CENTER)
      do k = kli, kui
        do j = jli, jui
          do i = ili, iui
            bsum = bsum + unk(ivar,i,j,k)**2
          enddo
        enddo
      enddo
    ! Release pointers:
    call Grid_releaseBlkPtr(lb,unk,CENTER)
    else
      do k = kli, kui
        do j = jli, jui
          do i = ili, iui
            bsum = bsum + work(i,j,k,lb,1)**2
          enddo
        enddo
      enddo
    endif
    lsum = lsum + bsum * cvol

  endif
enddo

call mpi_allreduce ( lsum, sum, 1, MPI_DOUBLE_PRECISION, & 
                     MPI_SUM, MPI_COMM_WORLD, ierr )
!call mpi_allreduce ( lvol, vol, 1, MPI_DOUBLE_PRECISION,
!                     MPI_SUM, MPI_COMM_WORLD, ierr )

norm = sqrt(sum)  ! definition in Briggs et al., _A Multigrid Tutorial_

!===============================================================================

return
end
