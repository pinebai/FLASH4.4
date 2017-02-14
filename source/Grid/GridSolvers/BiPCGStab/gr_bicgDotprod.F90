!!****if* source/Grid/GridSolvers/BiPCGStab/gr_bicgDotprod
!!
!! NAME
!!
!!  gr_bicgDotprod
!!
!! SYNOPSIS
!!
!!  call gr_bicgDotprod(integer(in) :: ivar1,
!!                      integer(in) :: ivar2,
!!                      real(inout) :: norm)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   ivar1 : 
!!
!!   ivar2 : 
!!
!!   norm : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

!*******************************************************************************

!  Routine:     bi_dotprod()

!  Description: Compute the do product of two variables. 

subroutine gr_bicgDotprod (ivar1,ivar2, norm)

!===============================================================================
#include "Flash.h"

  use gr_bicgData, ONLY: ili, jli, kli, iui, jui, kui

  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_getBlkPhysicalSize, Grid_getLocalNumBlks

  !temporary until we can think of better way
  !use dBaseDeclarations, ONLY: work
  use workspace, ONLY: work
  use tree, only : nodetype,newchild,lrefine

  use Timers_interface, ONLY: Timers_start, Timers_stop

  implicit none
#include "constants.h"
include "mpif.h"

  integer, intent(in) :: ivar1,ivar2
  real, intent(inout) :: norm

  integer :: lb, i, j, k, ierr, lnblocks
  real    :: lsum, bsum

  integer, parameter :: MAXDIMS = 3

  real, pointer, dimension(:,:,:,:), save :: unk

  integer, parameter :: nxb = NXB
  integer, parameter :: nyb = NYB
  integer, parameter :: nzb = NZB

!===============================================================================

!lvol = 0.
lsum = 0.

call Grid_getLocalNumBlks(lnblocks)

do lb = 1, lnblocks

  if (nodetype(lb) .eq. 1) then

    bsum = 0.

    if (ivar2 > 0) then
    ! Point to blocks center vars:
    call Grid_getBlkPtr(lb,unk,CENTER)
      do k = kli, kui
        do j = jli, jui
          do i = ili, iui
            bsum = bsum + unk(ivar1,i,j,k)*unk(ivar2,i,j,k)
          enddo
        enddo
      enddo
    ! Release pointers:
    call Grid_releaseBlkPtr(lb,unk,CENTER)

    else
    ! Point to blocks center vars:
    call Grid_getBlkPtr(lb,unk,CENTER)
      do k = kli, kui
        do j = jli, jui
          do i = ili, iui
            bsum = bsum + unk(ivar1,i,j,k)*work(i,j,k,lb,1)
          enddo
        enddo
      enddo
    ! Release pointers:
    call Grid_releaseBlkPtr(lb,unk,CENTER)
    endif

    lsum = lsum + bsum 

  endif
enddo

call Timers_start("mpi_allreduce")
call mpi_allreduce ( lsum, norm, 1, MPI_DOUBLE_PRECISION, & 
                     MPI_SUM, MPI_COMM_WORLD, ierr )
call Timers_stop("mpi_allreduce")

!===============================================================================

return
end subroutine gr_bicgDotprod
