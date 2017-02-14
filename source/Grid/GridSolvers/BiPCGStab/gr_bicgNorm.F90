!!****if* source/Grid/GridSolvers/BiPCGStab/gr_bicgNorm
!!
!! NAME
!!
!!  gr_bicgNorm
!!
!! SYNOPSIS
!!
!!  call gr_bicgNorm(integer(in) :: ivar,
!!                   real(inout) :: norm)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   ivar : 
!!
!!   norm : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

!*******************************************************************************

!  Routine:     bicg_norm()

!  Description: Compute the L2 norm of a given variable. 

subroutine gr_bicgNorm (ivar, norm)

!===============================================================================
#include "Flash.h"

  use gr_bicgData, ONLY: ili, jli, kli, iui, jui, kui

  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_getBlkPhysicalSize, Grid_getLocalNumBlks

  use workspace, ONLY: work
  use tree, only : nodetype,newchild,lrefine

  implicit none
#include "constants.h"
include "mpif.h"

  integer, intent(in) :: ivar
  real, intent(inout) :: norm

  integer :: lb, i, j, k, ierr, lnblocks
  real    :: lvol, vol, lsum, bsum, sum
  real    :: nbinv, cvol, bvol

  real, dimension(MDIM) :: size

  real, pointer, dimension(:,:,:,:), save :: unk

  integer, parameter :: nxb = NXB
  integer, parameter :: nyb = NYB
  integer, parameter :: nzb = NZB
  integer, parameter :: nguard = NGUARD
  integer, parameter :: ndim = NDIM
  integer, parameter :: k2d = K2D
  integer, parameter :: k3d = K3D

!===============================================================================

!!$lvol = 0.
lsum = 0.
nbinv = 1. / (real(nxb)*real(nyb)*real(nzb))

call Grid_getLocalNumBlks(lnblocks)

do lb = 1, lnblocks

  if (nodetype(lb) .eq. 1) then

    ! Get BlockSize:
    call Grid_getBlkPhysicalSize(lb,size)

    bvol = size(1)
    if (NDIM >= 2) bvol = bvol * size(2)
    if (NDIM == 3) bvol = bvol * size(3)
    cvol = bvol * nbinv
!!$    lvol = lvol + bvol
    bsum = 0.

    if (ivar > 0) then
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
end subroutine gr_bicgNorm
