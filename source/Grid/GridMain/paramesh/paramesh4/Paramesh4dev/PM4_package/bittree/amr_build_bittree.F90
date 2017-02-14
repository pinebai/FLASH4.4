!!****if* source/Grid/GridMain/paramesh/paramesh4/Paramesh4dev/PM4_package/bittree/amr_build_bittree
!!
!! NAME
!!
!!  amr_build_bittree
!!
!! SYNOPSIS
!!
!!  call amr_build_bittree()
!!
!! DESCRIPTION
!!
!!  Constructs the bittree from scratch by scanning all blocks in the mesh
!!  communicator.
!!
!!***
subroutine amr_build_bittree()
  use bittree, only : bittree_init,bittree_block_count,bittree_refine_init,bittree_refine_mark,&
       bittree_refine_apply,amr_verify_bittree 
  use paramesh_dimensions, only: ndim
  use Paramesh_comm_data, only: amr_mpi_meshComm
  use tree, only: lnblocks, coord, bsize, lrefine, nodetype, &
                  grid_xmin, grid_ymin, grid_zmin, &
                  grid_xmax, grid_ymax, grid_zmax
  use iso_c_binding, only: c_bool, c_int
  implicit none
#include "mpif.h"
  
  real :: gmin(ndim), gmax(ndim)
  integer :: top(3)
  integer(c_int) :: ijk(3)
  logical, allocatable :: topmask(:,:,:)
  integer(c_int) :: lev, b, pop0, pop1
  integer :: ierr
  
  gmin = reshape((/grid_xmin, grid_ymin, grid_zmin/), (/ndim/))
  gmax = reshape((/grid_xmax, grid_ymax, grid_zmax/), (/ndim/))
  ! all procs with blocks *should* agree on the nblockx/y/z values
  if(lnblocks > 0) then
    top(:) = 1
    top(1:ndim) = ishft(int(floor(0.5 + (gmax-gmin)/bsize(1:ndim,1))), 1-lrefine(1))
  else
    top(:) = 0
  end if
  call MPI_ALLREDUCE(MPI_IN_PLACE, top, 3, MPI_INTEGER, MPI_BOR, amr_mpi_meshComm, ierr)
  
  ! create topmask
  allocate(topmask(top(1),top(2),top(3)))
  topmask = .false.
  do b=1, lnblocks
    if(lrefine(b) == 1) then
      ijk(:) = 1
      ijk(1:ndim) = 1 + int((coord(1:ndim,b) - gmin)/bsize(1:ndim,b))
      topmask(ijk(1),ijk(2),ijk(3)) = .true.
    end if
  end do
  call MPI_ALLREDUCE(MPI_IN_PLACE, topmask, product(top), MPI_LOGICAL, MPI_LOR, amr_mpi_meshComm, ierr)
  
  ! create bittree
  call bittree_init(ndim, int(top,c_int), logical(topmask,c_bool))
  
  ! finer levels than top
  lev = 1
  pop0 = 0
  call bittree_block_count(pop1)
  do while(pop0 < pop1)
    call bittree_refine_init()
    do b=1, lnblocks
      if(lrefine(b) == lev .and. nodetype(b) > 1) then
        ijk(1:ndim) = int((coord(1:ndim,b) - gmin)/bsize(1:ndim,b), c_int)
        call bittree_refine_mark(lev-1, ijk)
      end if
    end do
    call bittree_refine_apply(amr_mpi_meshComm)
    pop0 = pop1
    call bittree_block_count(pop1)
    lev = lev + 1
  end do

  call amr_verify_bittree()
end subroutine
