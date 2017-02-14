!!****if* source/Grid/GridSolvers/Multigrid/gr_hgRecordNodeTypeState
!!
!! NAME
!!  gr_hgRecordNodeTypeState
!!
!! SYNOPSIS
!!  
!!  call gr_hgRecordNodeTypeState()
!!
!! DESCRIPTION
!!
!!  Save the original nodetypes using the gr_hgSaveNodetype and
!!  gr_hgSaveNewchild arrays, and record minimin and maximum
!!  levels of refinement of any existing blocks.
!!
!!
!!***


subroutine gr_hgRecordNodeTypeState

  use Logfile_interface, ONLY: Logfile_stamp

  use Grid_data, ONLY : gr_meshComm, gr_meshMe
  use gr_hgdata, ONLY: gr_hgSaveNodetype, gr_hgSaveNewchild, &
       gr_hgMeshRefineMin, gr_hgMeshRefineMax
  use tree, ONLY : lnblocks, newchild,nodetype, lrefine

  implicit none
#include "Flash_mpi.h"

  integer :: lb, ierr
  integer                            :: mylrefmin, mylrefmax

  do lb = 1, lnblocks
     gr_hgSaveNodetype(lb) = nodetype(lb)
     gr_hgSaveNewchild(lb) = newchild(lb)
  end do

  ! Determine minimum and maximum levels of refinement in the mesh.

  mylrefmin = HUGE(0)
  mylrefmax = 0

  do lb = 1, lnblocks
     mylrefmin = min(mylrefmin,lrefine(lb))
     mylrefmax = max(mylrefmax,lrefine(lb))
  end do
  
  call mpi_allreduce(mylrefmin, gr_hgMeshRefineMin, 1, MPI_INTEGER, MPI_MIN, &
       gr_meshComm, ierr)
  call mpi_allreduce(mylrefmax, gr_hgMeshRefineMax, 1, MPI_INTEGER, MPI_MAX, &
       gr_meshComm, ierr)
  
  !  mylrefmin = minval(lrefine(1:lnblocks))
  !  call mpi_allreduce(mylrefmin, gr_hgMeshRefineMin, 1, MPI_INTEGER, MPI_MIN, &
  !       gr_meshComm, ierr)
  
  !  mylrefmax = maxval(lrefine(1:lnblocks))
  !  call mpi_allreduce(mylrefmax, gr_hgMeshRefineMax, 1, MPI_INTEGER, MPI_MAX, &
  !       gr_meshComm, ierr)
  
  call Logfile_stamp(gr_hgMeshRefineMax,"[gr_hgRecordNodeTypeState] max refine level = ")

end subroutine gr_hgRecordNodeTypeState

