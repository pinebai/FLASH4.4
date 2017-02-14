!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/mpi_amr_boundary_block_info
!! NAME
!!
!!   mpi_amr_boundary_block_info
!!
!! SYNOPSIS
!!
!!   call mpi_amr_boundary_block_info (mype, nprocs)
!!
!!   call mpi_amr_boundary_block_info (integer, integer)
!!
!! ARGUMENTS
!!
!!   integer, intent(in) :: mype           
!!        The local processor number.
!!
!!   integer, intent(in) :: nprocs
!!        The number of MPI processes.
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!   mpif.h
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   mpi_morton
!!
!! CALLS
!!
!!   Does not call any other Paramesh routines.
!!
!! RETURNS
!!
!!   Does not return anything.  
!!
!! DESCRIPTION
!!
!!   This routine constructs a list of block faces which are external
!!   boundaries.
!!
!!
!! AUTHORS
!!
!!   Peter MacNeice (November 2002)
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      subroutine mpi_amr_boundary_block_info(mype,nprocs)

      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton
      use paramesh_comm_data

      implicit none

      include 'mpif.h'

      integer,intent(in) :: mype,nprocs

      integer,allocatable :: ib_global(:)
      integer,allocatable :: ib_count_send(:),ib_count_recv(:)
      integer,allocatable :: ib_global_recv(:)
      integer,allocatable :: ib_recv_index(:)
      integer,allocatable :: ib_send_index(:)

      integer :: ib,ib0,ib1,lb,ib_sum
      integer, allocatable :: ib_list(:,:)
      integer :: i,j,k,i1,i2,iproc,ierror

!-----------------------------------------------------------------

! Step 0.
! compute size for ib_list and then allocate it.
      ib0 = 0
      do lb = 1,lnblocks
        ib = minval(surr_blks(1,:,:,:,lb))
        if(ib.le.-20) then
          do k = 1,1+2*k3d
          do j = 1,1+2*k2d
          do i = 1,3
            ib1 = surr_blks(1,i,j,k,lb)
            if(ib1.le.-20) then
               ib0 = ib0 + 1
            end if
          enddo
          enddo
          enddo
        endif
      enddo
      if(allocated(ib_list)) deallocate(ib_list)
      allocate(ib_list(6,ib0))


! Step 1.
! construct on pe list of blocks next to boundary
      ib0 = 0
      do lb = 1,lnblocks
        ib = minval(surr_blks(1,:,:,:,lb))
        if(ib.le.-20) then
          do k = 1,1+2*k3d
          do j = 1,1+2*k2d
          do i = 1,3
            ib1 = surr_blks(1,i,j,k,lb)
            if(ib1.le.-20) then
            ib0 = ib0 + 1
            ib_list(1,ib0) = lb
            ib_list(2,ib0) = mype
            ib_list(3,ib0) = i
            ib_list(4,ib0) = j
            ib_list(5,ib0) = k
            ib_list(6,ib0) = ib
            endif
          enddo
          enddo
          enddo
        endif
      enddo

      if(.not.allocated(ib_global))  & 
     &                         allocate(ib_global(0:nprocs-1))
      if(.not.allocated(ib_global_recv))  & 
     &                         allocate(ib_global_recv(0:nprocs-1))
      if(.not.allocated(ib_count_send))  & 
     &                         allocate(ib_count_send(0:nprocs-1))
      if(.not.allocated(ib_count_recv))  & 
     &                         allocate(ib_count_recv(0:nprocs-1))
      if(.not.allocated(ib_recv_index))  & 
     &                  allocate(ib_recv_index(0:nprocs-1))
      if(.not.allocated(ib_send_index))  & 
     &                  allocate(ib_send_index(0:nprocs-1))
      ib_global_recv = 0
      ib_global = ib0

      call MPI_ALLTOALL (ib_global     ,1,MPI_INTEGER, & 
     &                   ib_global_recv,1,MPI_INTEGER, & 
     &                   amr_mpi_meshComm,ierror)
      ib_count_recv = ib_global_recv*6
      do iproc = 0,nprocs-1
        ib_count_send(iproc) = ib_count_recv(mype) 
      enddo

!
! Compute displacements in all to all message buffer, in bytes.
      ib_recv_index(0) = 0
      if(nprocs.gt.1) then
      do iproc = 1,nprocs-1
        ib_recv_index(iproc) = ib_recv_index(iproc-1)  & 
     &                        + ib_global_recv(iproc-1)*6
      enddo
      endif
      do iproc = 0,nprocs-1
        ib_send_index(iproc) = ib_recv_index(mype) 
      enddo



! Step 2.
! exchange no of boundary blocks on each processor
      call comm_int_sum_to_all(ib_sum,ib0)
      bc_block_neighs_length = ib_sum

! Step 3.
! allocate storage for global list
      if(allocated(bc_block_neighs)) deallocate(bc_block_neighs)
      allocate(bc_block_neighs(6,ib_sum))
      if(allocated(bc_block_neighs_send))  & 
     &                               deallocate(bc_block_neighs_send)
      allocate(bc_block_neighs_send(6,ib_sum))

! Step 4.
! exchange info between procs

! Put local data into is correct place on the list
      bc_block_neighs_send = -1
      bc_block_neighs = -5
      i1 = ib_recv_index(mype)/6 + 1
      i2 = i1 + ib0 - 1
      
      if(ib0.gt.0) bc_block_neighs_send(:,i1:i2) = ib_list(:,1:ib0)

! Exchange lists between all procs
      call MPI_ALLTOALLV (bc_block_neighs_send,ib_count_send, & 
     &                                    ib_send_index,MPI_INTEGER, & 
     &                    bc_block_neighs,ib_count_recv, & 
     &                                    ib_recv_index,MPI_INTEGER, & 
     &                    amr_mpi_meshComm,ierror)



      deallocate(ib_list)
      deallocate(ib_global)
      deallocate(ib_global_recv)
      deallocate(ib_count_send)
      deallocate(ib_count_recv)
      deallocate(ib_recv_index)
      deallocate(ib_send_index)

!
! Set status flag
      bc_block_neighs_status = 100

!-----------------------------------------------------------------

      return
      end subroutine mpi_amr_boundary_block_info
