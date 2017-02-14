!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f, unk_test, face[xyz]_test, edge[xyz]_test, unkn_test
#include "paramesh_preprocessor.fh"

!#define DEBUG

      subroutine amr_redist_blk(new_loc,nprocs,mype,lnblocks_old)

!!****f* mpi_source/amr_redist_blk
!! NAME
!!
!!   amr_redist_blk
!! 
!! SYNOPSIS
!!
!!   call amr_redist_blk (new_loc, nprocs, mype, lnblocks_old)
!!
!!   call amr_redist_blk (integer, integer, integer, integer, integer)
!!
!! ARGUMENTS      
!!
!!   integer, intent(inout) :: new_loc(:,:)
!!     Array which stores the new locations in memory where blocks are to be
!!     moved.  new_loc(1,:) indicates the on-processor location where the block
!!     is to reside and new_loc(2,:) indicates which processor the block data is
!!     to be moved to.
!!
!!   integer, intent(in) :: nprocs
!!     The number of processors being used.
!!   
!!   integer, intent (in) :: lnblocks_old
!!     The 'old' number of blocks in the calling processor before data redistribution.
!!
!!   integer, intent(in) :: mype     
!!      Current processor number
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   paramesh_interfaces
!!   paramesh_comm_data
!!
!! CALLS
!!
!!   fill_old_loc
!!   send_block_data
!!
!! RETURNS
!!
!!   Upon exit block data has been redistributed to processors according to
!!   a space filling morton curve.
!!
!! DESCRIPTION
!!
!!   This routine redistributes blocks to processors as part of the refinement/
!!   derefinement process.  This routine actually moves block data to their new
!!   locations if any refinements or derefinements have occured.
!!
!!   For each block a location in memory is passed in using the array new_loc.
!!   new_loc is the new location in memory where a block's data (e.g. in unk, 
!!   facevarx,y,z, ... etc.) is to be moved.  For instance, new_loc(1,lb) is the
!!   location in memory and new_loc(2,lb) is the processor to move the data to
!!   for block lb.
!!
!!   The routine is called by the subroutine amr_refine_derefine and should never 
!!   need to be called by a user's application.
!!
!! AUTHORS
!!
!!   Kevin Olson (1998)
!!   Bug fix contributed by Paul Ricker and Marcus Gross (2003)
!!
!!***

      use paramesh_dimensions
      use physicaldata
      use tree
      Use paramesh_comm_data

      use paramesh_interfaces, only : fill_old_loc

      implicit none

      include 'mpif.h'


      integer, intent(inout) :: new_loc(:,:)
      integer, intent(in)    :: nprocs,mype,lnblocks_old

      integer :: lb,ierr,errorcode
      logical :: free(maxblocks), moved(maxblocks), sent(maxblocks)
      logical :: repeat, repeatt
      integer :: old_loc(2,maxblocks_tr)
      integer :: nsend, nrecv
      integer :: reqr(maxblocks_tr)
      integer :: reqs(maxblocks_tr)
      integer :: statr(MPI_STATUS_SIZE,maxblocks_tr)
      integer :: stats(MPI_STATUS_SIZE,maxblocks_tr)
      integer :: nmoved, nit
      integer :: test(maxblocks), point_to(maxblocks)
      integer :: nm, nm2, nm2_old

      integer :: ireduce_datain(1),ireduce_dataout(1)
      logical :: lreduce_datain(1),lreduce_dataout(1)

      integer, save :: unk_int_type
      integer, allocatable :: unk_test(:,:,:,:)
      integer, save ::  is_unk,js_unk,ks_unk,ie_unk,je_unk,ke_unk

      integer, save :: facex_int_type, facey_int_type, facez_int_type
      integer, allocatable :: facex_test(:,:,:,:), facey_test(:,:,:,:),  &
     &                        facez_test(:,:,:,:)
      integer, save :: is_facex,js_facex,ks_facex
      integer, save :: ie_facex,je_facex,ke_facex
      integer, save :: is_facey,js_facey,ks_facey
      integer, save :: ie_facey,je_facey,ke_facey
      integer, save :: is_facez,js_facez,ks_facez
      integer, save :: ie_facez,je_facez,ke_facez

      integer, save :: edgex_int_type, edgey_int_type, edgez_int_type
      integer, allocatable :: edgex_test(:,:,:,:), edgey_test(:,:,:,:),  &
     &                        edgez_test(:,:,:,:)
      integer, save :: is_edgex,js_edgex,ks_edgex
      integer, save :: ie_edgex,je_edgex,ke_edgex
      integer, save :: is_edgey,js_edgey,ks_edgey
      integer, save :: ie_edgey,je_edgey,ke_edgey
      integer, save :: is_edgez,js_edgez,ks_edgez
      integer, save :: ie_edgez,je_edgez,ke_edgez

      integer, save :: unkn_int_type
      integer, allocatable :: unkn_test(:,:,:,:)
      integer, save ::  is_unkn,js_unkn,ks_unkn,ie_unkn,je_unkn,ke_unkn

      logical, save :: first = .true.
      integer :: type1, type2, type3
      integer :: udim(4), udim_tot(4), i
      integer :: nbytes

#ifdef DEBUG
      write(*,*) 'entering amr_redist_blk: pe ',mype
#endif /* DEBUG */

      if (first) then
      first = .false.

#ifdef REAL8
         nbytes = 8
#else
         nbytes = 4
#endif

      if (nvar > 0) then

      is_unk = nguard*npgs+1
      js_unk = nguard*k2d*npgs+1
      ks_unk = nguard*k3d*npgs+1
      ie_unk = nguard*npgs+nxb
      je_unk = nguard*k2d*npgs+nyb
      ke_unk = nguard*k3d*npgs+nzb

      ! The reordering script should take care of this
      allocate(unk_test(nvar,nxb,nyb,nzb))

      do i = 1,4
         udim_tot(i) = size(unk,dim=i) 
         udim(i) = size(unk_test,dim=i)
      end do

      deallocate(unk_test)

      ! DEFINE BLOCK INTERIOR
      call MPI_TYPE_VECTOR (udim(2),  & 
     &                      udim(1),  & 
     &                      udim_tot(1),  & 
     &                      amr_mpi_real,  & 
     &                      type1,  & 
     &                      ierr)
      call MPI_TYPE_HVECTOR (udim(3),  & 
     &                       1,  &
     &                       udim_tot(1)*udim_tot(2)*nbytes,  & 
     &                       type1,  & 
     &                       type2,  & 
     &                       ierr)
      call MPI_TYPE_HVECTOR (udim(4),  &
     &                       1,  &
     &                       udim_tot(1)*udim_tot(2)*udim_tot(3)*nbytes,  &
     &                       type2,  &
     &                       type3,  &
     &                       ierr)

      unk_int_type = type3

      call MPI_TYPE_COMMIT(unk_int_type,ierr)

      end if

      if (nfacevar > 0) then

      is_facex = nguard*npgs+1
      js_facex = nguard*k2d*npgs+1
      ks_facex = nguard*k3d*npgs+1
      ie_facex = nguard*npgs+nxb + 1
      je_facex = nguard*k2d*npgs+nyb
      ke_facex = nguard*k3d*npgs+nzb

      is_facey = nguard*npgs+1
      js_facey = nguard*k2d*npgs+1
      ks_facey = nguard*k3d*npgs+1
      ie_facey = nguard*npgs+nxb
      je_facey = nguard*k2d*npgs+nyb + k2d
      ke_facey = nguard*k3d*npgs+nzb

      is_facez = nguard*npgs+1
      js_facez = nguard*k2d*npgs+1
      ks_facez = nguard*k3d*npgs+1
      ie_facez = nguard*npgs+nxb 
      je_facez = nguard*k2d*npgs+nyb
      ke_facez = nguard*k3d*npgs+nzb + k3d

      ! The reordering script should take care of this
      allocate(facex_test(nfacevar,nxb+1,nyb,nzb))
      do i = 1,4
         udim_tot(i) = size(facevarx,dim=i) 
         udim(i) = size(facex_test,dim=i)
      end do
      deallocate(facex_test)

      ! DEFINE BLOCK INTERIOR
      call MPI_TYPE_VECTOR (udim(2),  & 
     &                      udim(1),  & 
     &                      udim_tot(1),  & 
     &                      amr_mpi_real,  & 
     &                      type1,  & 
     &                      ierr)
      call MPI_TYPE_HVECTOR (udim(3),  & 
     &                       1,  &
     &                       udim_tot(1)*udim_tot(2)*nbytes,  & 
     &                       type1,  & 
     &                       type2,  & 
     &                       ierr)
      call MPI_TYPE_HVECTOR (udim(4),  &
     &                       1,  &
     &                       udim_tot(1)*udim_tot(2)*udim_tot(3)*nbytes,  &
     &                       type2,  &
     &                       type3,  &
     &                       ierr)

      facex_int_type = type3

      call MPI_TYPE_COMMIT(facex_int_type,ierr)


      allocate(facey_test(nfacevar,nxb,nyb+k2d,nzb))
      do i = 1,4
         udim_tot(i) = size(facevary,dim=i) 
         udim(i) = size(facey_test,dim=i)
      end do
      deallocate(facey_test)

      ! DEFINE BLOCK INTERIOR
      call MPI_TYPE_VECTOR (udim(2),  & 
     &                      udim(1),  & 
     &                      udim_tot(1),  & 
     &                      amr_mpi_real,  & 
     &                      type1,  & 
     &                      ierr)
      call MPI_TYPE_HVECTOR (udim(3),  & 
     &                       1,  &
     &                       udim_tot(1)*udim_tot(2)*nbytes,  & 
     &                       type1,  & 
     &                       type2,  & 
     &                       ierr)
      call MPI_TYPE_HVECTOR (udim(4),  &
     &                       1,  &
     &                       udim_tot(1)*udim_tot(2)*udim_tot(3)*nbytes,  &
     &                       type2,  &
     &                       type3,  &
     &                       ierr)

      facey_int_type = type3

      call MPI_TYPE_COMMIT(facey_int_type,ierr)


      allocate(facez_test(nfacevar,nxb,nyb,nzb+k3d))
      do i = 1,4
         udim_tot(i) = size(facevarz,dim=i) 
         udim(i) = size(facez_test,dim=i)
      end do
      deallocate(facez_test)

      ! DEFINE BLOCK INTERIOR
      call MPI_TYPE_VECTOR (udim(2),  & 
     &                      udim(1),  & 
     &                      udim_tot(1),  & 
     &                      amr_mpi_real,  & 
     &                      type1,  & 
     &                      ierr)
      call MPI_TYPE_HVECTOR (udim(3),  & 
     &                       1,  &
     &                       udim_tot(1)*udim_tot(2)*nbytes,  & 
     &                       type1,  & 
     &                       type2,  & 
     &                       ierr)
      call MPI_TYPE_HVECTOR (udim(4),  &
     &                       1,  &
     &                       udim_tot(1)*udim_tot(2)*udim_tot(3)*nbytes,  &
     &                       type2,  &
     &                       type3,  &
     &                       ierr)

      facez_int_type = type3

      call MPI_TYPE_COMMIT(facez_int_type,ierr)

      end if

      if (nvaredge > 0) then

      is_edgex = nguard*npgs+1
      js_edgex = nguard*k2d*npgs+1
      ks_edgex = nguard*k3d*npgs+1
      ie_edgex = nguard*npgs+nxb
      je_edgex = nguard*k2d*npgs+nyb + k2d
      ke_edgex = nguard*k3d*npgs+nzb + k3d

      is_edgey = nguard*npgs+1
      js_edgey = nguard*k2d*npgs+1
      ks_edgey = nguard*k3d*npgs+1
      ie_edgey = nguard*npgs+nxb + 1
      je_edgey = nguard*k2d*npgs+nyb
      ke_edgey = nguard*k3d*npgs+nzb + k3d

      is_edgez = nguard*npgs+1
      js_edgez = nguard*k2d*npgs+1
      ks_edgez = nguard*k3d*npgs+1
      ie_edgez = nguard*npgs+nxb + 1
      je_edgez = nguard*k2d*npgs+nyb + k2d
      ke_edgez = nguard*k3d*npgs+nzb

      ! The reordering script should take care of this
      allocate(edgex_test(nvaredge,nxb,nyb+k2d,nzb+k3d))
      do i = 1,4
         udim_tot(i) = size(unk_e_x,dim=i) 
         udim(i) = size(edgex_test,dim=i)
      end do
      deallocate(edgex_test)

      ! DEFINE BLOCK INTERIOR
      call MPI_TYPE_VECTOR (udim(2),  & 
     &                      udim(1),  & 
     &                      udim_tot(1),  & 
     &                      amr_mpi_real,  & 
     &                      type1,  & 
     &                      ierr)
      call MPI_TYPE_HVECTOR (udim(3),  & 
     &                       1,  &
     &                       udim_tot(1)*udim_tot(2)*nbytes,  & 
     &                       type1,  & 
     &                       type2,  & 
     &                       ierr)
      call MPI_TYPE_HVECTOR (udim(4),  &
     &                       1,  &
     &                       udim_tot(1)*udim_tot(2)*udim_tot(3)*nbytes,  &
     &                       type2,  &
     &                       type3,  &
     &                       ierr)

      edgex_int_type = type3

      call MPI_TYPE_COMMIT(edgex_int_type,ierr)


      allocate(edgey_test(nvaredge,nxb+1,nyb,nzb+k3d))
      do i = 1,4
         udim_tot(i) = size(unk_e_y,dim=i) 
         udim(i) = size(edgey_test,dim=i)
      end do
      deallocate(edgey_test)

      ! DEFINE BLOCK INTERIOR
      call MPI_TYPE_VECTOR (udim(2),  & 
     &                      udim(1),  & 
     &                      udim_tot(1),  & 
     &                      amr_mpi_real,  & 
     &                      type1,  & 
     &                      ierr)
      call MPI_TYPE_HVECTOR (udim(3),  & 
     &                       1,  &
     &                       udim_tot(1)*udim_tot(2)*nbytes,  & 
     &                       type1,  & 
     &                       type2,  & 
     &                       ierr)
      call MPI_TYPE_HVECTOR (udim(4),  &
     &                       1,  &
     &                       udim_tot(1)*udim_tot(2)*udim_tot(3)*nbytes,  &
     &                       type2,  &
     &                       type3,  &
     &                       ierr)

      edgey_int_type = type3

      call MPI_TYPE_COMMIT(edgey_int_type,ierr)


      allocate(edgez_test(nvaredge,nxb+1,nyb+k2d,nzb))
      do i = 1,4
         udim_tot(i) = size(unk_e_z,dim=i) 
         udim(i) = size(edgez_test,dim=i)
      end do
      deallocate(edgez_test)

      ! DEFINE BLOCK INTERIOR
      call MPI_TYPE_VECTOR (udim(2),  & 
     &                      udim(1),  & 
     &                      udim_tot(1),  & 
     &                      amr_mpi_real,  & 
     &                      type1,  & 
     &                      ierr)
      call MPI_TYPE_HVECTOR (udim(3),  & 
     &                       1,  &
     &                       udim_tot(1)*udim_tot(2)*nbytes,  & 
     &                       type1,  & 
     &                       type2,  & 
     &                       ierr)
      call MPI_TYPE_HVECTOR (udim(4),  &
     &                       1,  &
     &                       udim_tot(1)*udim_tot(2)*udim_tot(3)*nbytes,  &
     &                       type2,  &
     &                       type3,  &
     &                       ierr)

      edgez_int_type = type3

      call MPI_TYPE_COMMIT(edgez_int_type,ierr)

      end if

      if (nvarcorn > 0) then

      is_unkn = nguard*npgs+1
      js_unkn = nguard*k2d*npgs+1
      ks_unkn = nguard*k3d*npgs+1
      ie_unkn = nguard*npgs+nxb + 1
      je_unkn = nguard*k2d*npgs+nyb + k2d
      ke_unkn = nguard*k3d*npgs+nzb + k3d

      ! The reordering script should take care of this
      allocate(unkn_test(nvarcorn,nxb+1,nyb+k2d,nzb+k3d))

      do i = 1,4
         udim_tot(i) = size(unk_n,dim=i) 
         udim(i) = size(unkn_test,dim=i)
      end do
      deallocate(unkn_test)

      ! DEFINE BLOCK INTERIOR
      call MPI_TYPE_VECTOR (udim(2),  & 
     &                      udim(1),  & 
     &                      udim_tot(1),  & 
     &                      amr_mpi_real,  & 
     &                      type1,  & 
     &                      ierr)
      call MPI_TYPE_HVECTOR (udim(3),  & 
     &                       1,  &
     &                       udim_tot(1)*udim_tot(2)*nbytes,  & 
     &                       type1,  & 
     &                       type2,  & 
     &                       ierr)
      call MPI_TYPE_HVECTOR (udim(4),  &
     &                       1,  &
     &                       udim_tot(1)*udim_tot(2)*udim_tot(3)*nbytes,  &
     &                       type2,  &
     &                       type3,  &
     &                       ierr)

      unkn_int_type = type3

      call MPI_TYPE_COMMIT(unkn_int_type,ierr)

      end if

      end if


! 1) compute old_loc
      call fill_old_loc (new_loc,old_loc,nprocs,mype)

      nrecv = 0
      nsend = 0


!--------------
! treat unk

      if(nvar.gt.0) then

! Post all receives for unk
         do lb = 1,new_lnblocks
            if (.not.newchild(lb)) then
               if (old_loc(2,lb).ne.mype) then
                  nrecv = nrecv + 1
                  call MPI_IRECV (unk(1,is_unk,js_unk,ks_unk,lb), & 
     &                            1, & 
     &                            unk_int_type, & 
     &                            old_loc(2,lb), & 
     &                            lb, & 
     &                            amr_mpi_meshComm, & 
     &                            reqr(nrecv), & 
     &                            ierr)

               end if
            end if
         end do

      end if

!--------------------
! Treat Facevariables
!--------------------

      if (nfacevar.gt.0) then

!---------------
! Treat facevarx
!---------------

         do lb = 1,new_lnblocks
            if (.not.newchild(lb)) then
               if (old_loc(2,lb).ne.mype) then
                  nrecv = nrecv + 1
                  call MPI_IRECV (facevarx(1,is_facex,js_facex,ks_facex,lb), & 
     &                            1, & 
     &                            facex_int_type, & 
     &                            old_loc(2,lb), & 
     &                            lb+2*maxblocks, & 
     &                            amr_mpi_meshComm, & 
     &                            reqr(nrecv), & 
     &                            ierr)
               end if
            end if
         end do

!---------------
! Treat facevary
!---------------

         if (ndim >= 2) then
         do lb = 1,new_lnblocks
            if (.not.newchild(lb)) then
               if (old_loc(2,lb).ne.mype) then
                  nrecv = nrecv + 1
                  call MPI_IRECV (facevary(1,is_facey,js_facey,ks_facey,lb), & 
     &                            1, & 
     &                            facey_int_type, & 
     &                            old_loc(2,lb), & 
     &                            lb+3*maxblocks, & 
     &                            amr_mpi_meshComm, & 
     &                            reqr(nrecv), & 
     &                            ierr)
               end if
            end if
         end do
         end if

!---------------
! Treat Facevarz
!---------------

         if (ndim == 3) then
         do lb = 1,new_lnblocks
            if (.not.newchild(lb)) then
               if (old_loc(2,lb).ne.mype) then
                  nrecv = nrecv + 1
                  call MPI_IRECV (facevarz(1,is_facez,js_facez,ks_facez,lb), & 
     &                            1, & 
     &                            facez_int_type, & 
     &                            old_loc(2,lb), & 
     &                            lb+4*maxblocks, & 
     &                            amr_mpi_meshComm, & 
     &                            reqr(nrecv), & 
     &                            ierr)
               end if
            end if
         end do
         end if

      end if

!--------------------
! Treat Edge variables
!--------------------

      if (nvaredge.gt.0) then

!---------------
! Treat unk_e_x
!---------------

         do lb = 1,new_lnblocks
            if (.not.newchild(lb)) then
               if (old_loc(2,lb).ne.mype) then
                  nrecv = nrecv + 1
                  call MPI_IRECV (unk_e_x(1,is_edgex,js_edgex,ks_edgex,lb), & 
     &                            1, & 
     &                            edgex_int_type, & 
     &                            old_loc(2,lb), & 
     &                            lb+5*maxblocks, & 
     &                            amr_mpi_meshComm, & 
     &                            reqr(nrecv), & 
     &                            ierr)
               end if
            end if
         end do

!---------------
! Treat unk_e_y
!---------------

         if (ndim >= 2) then
         do lb = 1,new_lnblocks
            if (.not.newchild(lb)) then
               if (old_loc(2,lb).ne.mype) then
                  nrecv = nrecv + 1
                  call MPI_IRECV (unk_e_y(1,is_edgey,js_edgey,ks_edgey,lb), & 
     &                            1, & 
     &                            edgey_int_type, & 
     &                            old_loc(2,lb), & 
     &                            lb+6*maxblocks, & 
     &                            amr_mpi_meshComm, & 
     &                            reqr(nrecv), & 
     &                            ierr)
               end if
            end if
         end do
         end if
!---------------
! Treat unk_e_z
!---------------

         if (ndim == 3) then
         do lb = 1,new_lnblocks
            if (.not.newchild(lb)) then
               if (old_loc(2,lb).ne.mype) then
                  nrecv = nrecv + 1
                  call MPI_IRECV (unk_e_z(1,is_edgez,js_edgez,ks_edgez,lb), & 
     &                            1, & 
     &                            edgez_int_type, & 
     &                            old_loc(2,lb), & 
     &                            lb+7*maxblocks, & 
     &                            amr_mpi_meshComm, & 
     &                            reqr(nrecv), & 
     &                            ierr)
               end if
            end if
         end do
         end if

      endif

!--------------
! treat unk_n

      if(nvarcorn.gt.0) then

! Post all receives for unk_n
         do lb = 1,new_lnblocks
            if (.not.newchild(lb)) then
               if (old_loc(2,lb).ne.mype) then
                  nrecv = nrecv + 1
                  call MPI_IRECV (unk_n(1,is_unkn,js_unkn,ks_unkn,lb), & 
     &                            1, & 
     &                            unkn_int_type, & 
     &                            old_loc(2,lb), & 
     &                            lb+8*maxblocks, & 
     &                            amr_mpi_meshComm, & 
     &                            reqr(nrecv), & 
     &                            ierr)
               end if
            end if
         end do

      end if

!--------------

      
      moved(:) = .false.
      moved(lnblocks_old+1:maxblocks) = .true.
      free(:) = .false.
      free(lnblocks_old+1:maxblocks) = .true.
      sent(:) = .false.
      repeat = .TRUE.
      nmoved = 0 
      test(:) = 0
      point_to(:) = 0
      
      nit = 0
      nm2 = 0
      nm2_old = 1
      do while (repeat.and.nit<=100) 
         
         do lb = 1,max(lnblocks_old,new_lnblocks)
            call send_block_data (lb, new_loc, old_loc, free,  & 
     &                            moved, sent, & 
     &                            lnblocks_old, mype, nmoved, & 
     &                            test, point_to, & 
     &                            reqs, nsend, unk_int_type,  &
     &                            facex_int_type, facey_int_type, &
     &                            facez_int_type, edgex_int_type, &
     &                            edgey_int_type, edgez_int_type, &
     &                            unkn_int_type)
         end do
         repeat = any(.not.moved(:))
         lreduce_datain(1) = repeat
         call mpi_logical_allreduce( & 
     &           lreduce_datain(1),lreduce_dataout(1), & 
     &           1,MPI_LOGICAL, & 
     &           MPI_LOR,amr_mpi_meshComm,ierr)
         repeatt = lreduce_dataout(1)
         repeat = repeatt
         
         nm2_old = nm2
         nm = count(.not.moved(:))
         ireduce_datain(1) = nm
         call mpi_int_allreduce( & 
     &        ireduce_datain(1),ireduce_dataout(1), & 
     &        1,MPI_INTEGER, & 
     &        MPI_SUM,amr_mpi_meshComm,ierr)
         nm2 = ireduce_dataout(1)
         if (mype.eq.0) then
            print *,' iteration, no. not moved = ',nit,nm2
         end if
         
         nit = nit + 1
         
      end do
      
      if (nm2_old.eq.nm2.and.nm2.ne.0.and.nit>=100) then
         if (mype.eq.0) then
          print *,' ERROR: could not move all blocks in amr_redist_blk'
          print *,' Try increasing maxblocks or use more processors'
          print *,' nm2_old, nm2 = ',nm2_old,nm2
          print *,' ABORTING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         end if
         call MPI_ABORT(amr_mpi_meshComm,errorcode,ierr)
      end if
      
      if (nrecv.gt.0) then
         call MPI_WAITALL(nrecv,reqr,statr,ierr)
      end if

! NOTE: Bug fix by Paul Ricker and Marcus Gross (5/2003).  Added Waitall 
! to isends so that SGI MPI buffers do not overflow.

      if (nsend.gt.0) then
         call MPI_WAITALL(nsend,reqs,stats,ierr)
      end if

#ifdef DEBUG
      write(*,*) 'exiting amr_redist_blk: pe ',mype
#endif /* DEBUG */
      
      return
      end subroutine amr_redist_blk
      
!----------------------------------------------------------------------------
      
      subroutine send_block_data (lb, new_loc, old_loc, free, & 
     &                            moved, sent, & 
     &                            lnblocks_old, mype, nmoved, & 
     &                            test, point_to, & 
     &                            reqs, nsend, unk_int_type,  &
     &                            facex_int_type, facey_int_type, &
     &                            facez_int_type, edgex_int_type, &
     &                            edgey_int_type, edgez_int_type, &
     &                            unkn_int_type)
         
!!****f* mpi_source/send_block_data
!! NAME
!!
!!   send_block_data
!! 
!! SYNOPSIS
!!
!!   call send_block_data (lb, new_loc, old_loc, free, moved, sent, 
!!                         lnblocks_old, mype, nmoved, test, point_to,
!!                         reqs, nsend, unk_int_type,
!!                         facex_int_type, facey_int_type,
!!                         facez_int_type, edgex_int_type,
!!                         edgey_int_type, edgez_int_type,
!!                         unkn_int_type)
!!
!!   call send_block_data (integer, integer, integer, logical, logical, logical,
!!                         integer, integer, integer, integer, integer,
!!                         integer, integer, integer,
!!                         integer, integer, integer)
!!
!! ARGUMENTS      
!!
!!   integer :: lb
!!     Block being sent
!!
!!   integer :: new_loc(2,maxblocks_tr)
!!     Array which stores the new locations in memory where blocks are to be
!!     moved.  new_loc(1,:) indicates the on-processor location where the block
!!     is to reside and new_loc(2,:) indicates which processor the block data is
!!     to be moved to.
!!
!!   integer :: old_loc(2,maxblocks_tr)
!!     Similar to new_loc above except stores the location in memory where a 
!!     block orginates from.  Necessary for hand shaking during message passing.
!!
!!   logical :: free(maxblocks)
!!     A logical placeholder which indicates if a position in the list of blocks
!!     is free to receive data.
!!
!!   logical :: moved(maxblocks)
!!     A logical which indicates if data has moved from a location in the list of 
!!     blocks.
!!
!!   logical :: sent(maxblocks)
!!     A logical which indicates if data has been sent from a location in the list 
!!     of blocks.
!!
!!   integer :: lnblocks_old
!!     The 'old' number of blocks in the calling processor before data redistribution.
!!
!!   integer :: mype     
!!     Current processor number
!!
!!   integer :: nmoved
!!     Indicates the number of blocks which have been moved during the algorithm. 
!!
!!   integer :: test(maxblocks)
!!     ???
!!
!!   integer :: point_to(maxblocks)
!!     ???
!!   
!!   integer :: reqs(maxblocks_tr)
!!     List of MPI requests.
!!
!!   integer :: nsend
!!     Number MPI sends issued.
!!
!!   integer :: unk_int_type
!!     An MPI derived type defining the interior region of block not
!!     including guardcells (unk).
!!
!!   integer :: facex_int_type, facey_int_type, facez_int_type
!!     An MPI derived type defining the interior region of block not
!!     including guardcells (face variables)
!!
!!   integer :: edgex_int_type, edgey_int_type, edgez_int_type
!!     An MPI derived type defining the interior region of block not
!!     including guardcells (edge variables)
!!
!!   integer :: unkn_int_type
!!     An MPI derived type defining the interior region of block not
!!     including guardcells (unk_n).
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!
!! CALLS
!!
!!
!! RETURNS
!!
!!
!! DESCRIPTION
!!
!!   This routine is part of the redistribution of block data and manages the 
!!   sending of block data to their new locations.  The block data is sent to
!!   its final new location (new_loc) if that location is unoccupied by another
!!   block's data, it is sent to another free location, or nothing is sent.
!!
!!   The routine is only called by the subroutine amr_redist_blk and should never
!!   need to be called by a user's application.
!!
!! AUTHORS
!!
!!   Kevin Olson (1998)
!!
!!***

      use paramesh_dimensions
      use physicaldata
      use tree
      use paramesh_comm_data

      implicit none

      include 'mpif.h'


      integer :: new_loc(2,maxblocks_tr), old_loc(2,maxblocks_tr)
      logical :: free(maxblocks), moved(maxblocks), sent(maxblocks)
      integer :: lb, lnblocks_old, mype
      logical :: success
      integer :: status(MPI_STATUS_SIZE)
      integer :: reqs(maxblocks_tr), nsend
      integer :: ierr, nmoved, lb2
      integer :: point_to(maxblocks),test(maxblocks)

      integer :: unk_int_type
      integer :: facex_int_type, facey_int_type, facez_int_type
      integer :: edgex_int_type, edgey_int_type, edgez_int_type
      integer :: unkn_int_type

      logical, save :: first = .true.

      integer,save ::  is_unk,js_unk,ks_unk,ie_unk,je_unk,ke_unk

      integer,save ::  is_facex,js_facex,ks_facex,ie_facex,je_facex,ke_facex
      integer,save ::  is_facey,js_facey,ks_facey,ie_facey,je_facey,ke_facey
      integer,save ::  is_facez,js_facez,ks_facez,ie_facez,je_facez,ke_facez

      integer,save ::  is_edgex,js_edgex,ks_edgex,ie_edgex,je_edgex,ke_edgex
      integer,save ::  is_edgey,js_edgey,ks_edgey,ie_edgey,je_edgey,ke_edgey
      integer,save ::  is_edgez,js_edgez,ks_edgez,ie_edgez,je_edgez,ke_edgez

      integer,save ::  is_unkn,js_unkn,ks_unkn,ie_unkn,je_unkn,ke_unkn

      if (first) then

      first = .false.

      is_unk = nguard*npgs+1
      js_unk = nguard*k2d*npgs+1
      ks_unk = nguard*k3d*npgs+1
      ie_unk = nguard*npgs+nxb
      je_unk = nguard*k2d*npgs+nyb
      ke_unk = nguard*k3d*npgs+nzb

      is_facex = nguard*npgs+1
      js_facex = nguard*k2d*npgs+1
      ks_facex = nguard*k3d*npgs+1
      ie_facex = nguard*npgs+nxb + 1
      je_facex = nguard*k2d*npgs+nyb
      ke_facex = nguard*k3d*npgs+nzb

      is_facey = nguard*npgs+1
      js_facey = nguard*k2d*npgs+1
      ks_facey = nguard*k3d*npgs+1
      ie_facey = nguard*npgs+nxb 
      je_facey = nguard*k2d*npgs+nyb + k2d
      ke_facey = nguard*k3d*npgs+nzb

      is_facez = nguard*npgs+1
      js_facez = nguard*k2d*npgs+1
      ks_facez = nguard*k3d*npgs+1
      ie_facez = nguard*npgs+nxb 
      je_facez = nguard*k2d*npgs+nyb
      ke_facez = nguard*k3d*npgs+nzb + k3d

      is_edgex = nguard*npgs+1
      js_edgex = nguard*k2d*npgs+1
      ks_edgex = nguard*k3d*npgs+1
      ie_edgex = nguard*npgs+nxb
      je_edgex = nguard*k2d*npgs+nyb + k2d
      ke_edgex = nguard*k3d*npgs+nzb + k3d

      is_edgey = nguard*npgs+1
      js_edgey = nguard*k2d*npgs+1
      ks_edgey = nguard*k3d*npgs+1
      ie_edgey = nguard*npgs+nxb + 1
      je_edgey = nguard*k2d*npgs+nyb
      ke_edgey = nguard*k3d*npgs+nzb + k3d

      is_edgez = nguard*npgs+1
      js_edgez = nguard*k2d*npgs+1
      ks_edgez = nguard*k3d*npgs+1
      ie_edgez = nguard*npgs+nxb + 1
      je_edgez = nguard*k2d*npgs+nyb + k2d
      ke_edgez = nguard*k3d*npgs+nzb

      is_unkn = nguard*npgs+1
      js_unkn = nguard*k2d*npgs+1
      ks_unkn = nguard*k3d*npgs+1
      ie_unkn = nguard*npgs+nxb + 1
      je_unkn = nguard*k2d*npgs+nyb + k2d
      ke_unkn = nguard*k3d*npgs+nzb + k3d

      end if

      if (new_loc(1,lb).eq.lb.and.new_loc(2,lb).eq.mype) then
         if (.not.moved(lb)) moved(lb) = .true.
         return
      end if

      if (lb.le.max(lnblocks_old,new_lnblocks)) then

         if (lb.le.lnblocks_old) then
           if (new_loc(2,lb).ne.mype) then
            success = .false.
            call MPI_IPROBE (new_loc(2,lb), & 
     &                       maxblocks+new_loc(1,lb), & 
     &                       amr_mpi_meshComm, & 
     &                       success, & 
     &                       status, & 
     &                       ierr)
            if (.not.moved(lb).and.success) then
               call MPI_RECV (success, & 
     &                        1, & 
     &                        MPI_LOGICAL, & 
     &                        new_loc(2,lb), & 
     &                        maxblocks+new_loc(1,lb), & 
     &                        amr_mpi_meshComm, & 
     &                        status, & 
     &                        ierr)
               if (free(lb)) then
                if (nvar.gt.0) then
                call MPI_SSEND ( & 
     &   unk(1,is_unk,js_unk,ks_unk,point_to(lb)), & 
     &                          1, & 
     &                          unk_int_type, & 
     &                          new_loc(2,lb), & 
     &                          new_loc(1,lb), & 
     &                          amr_mpi_meshComm, & 
     &                          ierr)
                endif
! send facevariables
                if (nfacevar.gt.0) then
                call MPI_SSEND (facevarx(1,is_facex,js_facex,ks_facex, &
     &                          point_to(lb)), & 
     &                          1, & 
     &                          facex_int_type, & 
     &                          new_loc(2,lb), & 
     &                          new_loc(1,lb)+2*maxblocks, & 
     &                          amr_mpi_meshComm, & 
     &                          ierr)
                if (ndim >= 2) then
                call MPI_SSEND (facevary(1,is_facey,js_facey,ks_facey,  &
     &                          point_to(lb)), & 
     &                          1, & 
     &                          facey_int_type, & 
     &                          new_loc(2,lb), & 
     &                          new_loc(1,lb)+3*maxblocks, & 
     &                          amr_mpi_meshComm, & 
     &                          ierr)
                end if
                if (ndim == 3) then
                call MPI_SSEND (facevarz(1,is_facez,js_facez,ks_facez,  &
     &                          point_to(lb)), & 
     &                          1, & 
     &                          facez_int_type, & 
     &                          new_loc(2,lb), & 
     &                          new_loc(1,lb)+4*maxblocks, & 
     &                          amr_mpi_meshComm, & 
     &                          ierr)
                end if

                end if

! send edge variables
                if (nvaredge.gt.0) then
                call MPI_SSEND (unk_e_x(1,is_edgex,js_edgex,ks_edgex,  &
     &                          point_to(lb)), & 
     &                          1, & 
     &                          edgex_int_type, & 
     &                          new_loc(2,lb), & 
     &                          new_loc(1,lb)+5*maxblocks, & 
     &                          amr_mpi_meshComm, & 
     &                          ierr)
                if (ndim >= 2) then
                call MPI_SSEND (unk_e_y(1,is_edgey,js_edgey,ks_edgey, &
     &                          point_to(lb)), & 
     &                          1, & 
     &                          edgey_int_type, & 
     &                          new_loc(2,lb), & 
     &                          new_loc(1,lb)+6*maxblocks, & 
     &                          amr_mpi_meshComm, & 
     &                          ierr)
                end if
                if (ndim == 3) then
                call MPI_SSEND (unk_e_z(1,is_edgez,js_edgez,ks_edgez,  &
     &                          point_to(lb)), & 
     &                          1, & 
     &                          edgez_int_type, & 
     &                          new_loc(2,lb), & 
     &                          new_loc(1,lb)+7*maxblocks, & 
     &                          amr_mpi_meshComm, & 
     &                          ierr)
                end if
                end if

! send corner variables
                if (nvarcorn.gt.0) then
                call MPI_SSEND (unk_n(1,is_unkn,js_unkn,ks_unkn,point_to(lb)), & 
     &                          1, & 
     &                          unkn_int_type, & 
     &                          new_loc(2,lb), & 
     &                          new_loc(1,lb)+8*maxblocks, & 
     &                          amr_mpi_meshComm, & 
     &                          ierr)
                end if

                test(point_to(lb)) = -1
               else
                if (nvar.gt.0) then
                call MPI_SSEND (unk(1,is_unk,js_unk,ks_unk,lb), & 
     &                          1, & 
     &                          unk_int_type, & 
     &                          new_loc(2,lb), & 
     &                          new_loc(1,lb), & 
     &                          amr_mpi_meshComm, & 
     &                          ierr)
                end if
! send facevariables
                if (nfacevar.gt.0) then
                call MPI_SSEND (facevarx(1,is_facex,js_facex,ks_facex,lb), & 
     &                          1, & 
     &                          facex_int_type, & 
     &                          new_loc(2,lb), & 
     &                          new_loc(1,lb)+2*maxblocks, & 
     &                          amr_mpi_meshComm, & 
     &                          ierr)
                if (ndim >= 2) then
                call MPI_SSEND (facevary(1,is_facey,js_facey,ks_facey,lb), & 
     &                          1, & 
     &                          facey_int_type, & 
     &                          new_loc(2,lb), & 
     &                          new_loc(1,lb)+3*maxblocks, & 
     &                          amr_mpi_meshComm, & 
     &                          ierr)
                end if
                if (ndim == 3) then
                call MPI_SSEND (facevarz(1,is_facez,js_facez,ks_facez,lb), & 
     &                          1, & 
     &                          facez_int_type, & 
     &                          new_loc(2,lb), & 
     &                          new_loc(1,lb)+4*maxblocks, & 
     &                          amr_mpi_meshComm, & 
     &                          ierr)
                end if
                end if

! send edge variables
                if (nvaredge.gt.0) then
                call MPI_SSEND (unk_e_x(1,is_edgex,js_edgex,ks_edgex,lb), & 
     &                          1, & 
     &                          edgex_int_type, & 
     &                          new_loc(2,lb), & 
     &                          new_loc(1,lb)+5*maxblocks, & 
     &                          amr_mpi_meshComm, & 
     &                          ierr)
                if (ndim >= 2) then
                call MPI_SSEND (unk_e_y(1,is_edgey,js_edgey,ks_edgey,lb), & 
     &                          1, & 
     &                          edgey_int_type, & 
     &                          new_loc(2,lb), & 
     &                          new_loc(1,lb)+6*maxblocks, & 
     &                          amr_mpi_meshComm, & 
     &                          ierr)
                end if
                if (ndim == 3) then
                call MPI_SSEND (unk_e_z(1,is_edgez,js_edgez,ks_edgez,lb), & 
     &                          1, & 
     &                          edgez_int_type, & 
     &                          new_loc(2,lb), & 
     &                          new_loc(1,lb)+7*maxblocks, & 
     &                          amr_mpi_meshComm, & 
     &                          ierr)
                end if
                end if

! send corner variables
                if (nvarcorn.gt.0) then
                call MPI_SSEND (unk_n(1,is_unkn,js_unkn,ks_unkn,lb), & 
     &                          1, & 
     &                          unkn_int_type, & 
     &                          new_loc(2,lb), & 
     &                          new_loc(1,lb)+8*maxblocks, & 
     &                          amr_mpi_meshComm, & 
     &                          ierr)
                end if

                free(lb) = .true.
               end if
               moved(lb) = .true.
            end if
           else
            if (.not.moved(lb).and.free(new_loc(1,lb))) then
             if (free(lb)) then

               if (nvar.gt.0) then
         unk(:,is_unk:ie_unk,js_unk:je_unk,ks_unk:ke_unk,new_loc(1,lb)) =  & 
     &    unk(:,is_unk:ie_unk,js_unk:je_unk,ks_unk:ke_unk,point_to(lb))
               end if

! move facevars
               if (nfacevar.gt.0) then
                  facevarx(:,is_facex:ie_facex,js_facex:je_facex,  &
     &                       ks_facex:ke_facex,new_loc(1,lb)) =  & 
     &                 facevarx(:,is_facex:ie_facex,js_facex:je_facex, &
     &                            ks_facex:ke_facex,point_to(lb))
                  if (ndim >= 2) then
                  facevary(:,is_facey:ie_facey,js_facey:je_facey, &
     &                       ks_facey:ke_facey,new_loc(1,lb)) =  & 
     &                 facevary(:,is_facey:ie_facey,js_facey:je_facey,  &
     &                            ks_facey:ke_facey,point_to(lb))
                  end if
                  if (ndim == 3) then
                  facevarz(:,is_facez:ie_facez,js_facez:je_facez,  &
     &                       ks_facez:ke_facez,new_loc(1,lb)) =  & 
     &                 facevarz(:,is_facez:ie_facez,js_facez:je_facez,  &
     &                            ks_facez:ke_facez,point_to(lb))
                  end if
               end if

! move edgevars
               if (nvaredge.gt.0) then
                  unk_e_x(:,is_edgex:ie_edgex,js_edgex:je_edgex,  &
     &                      ks_edgex:ke_edgex,new_loc(1,lb)) = & 
     &                 unk_e_x(:,is_edgex:ie_edgex,js_edgex:je_edgex,  &
     &                           ks_edgex:ke_edgex,point_to(lb))
                  if (ndim >= 2) then
                  unk_e_y(:,is_edgey:ie_edgey,js_edgey:je_edgey,  &
     &                      ks_edgey:ke_edgey,new_loc(1,lb)) = & 
     &                 unk_e_y(:,is_edgey:ie_edgey,js_edgey:je_edgey,  &
     &                           ks_edgey:ke_edgey,point_to(lb))
                  end if
                  if (ndim == 3) then
                  unk_e_z(:,is_edgez:ie_edgez,js_edgez:je_edgez,  &
     &                      ks_edgez:ke_edgez,new_loc(1,lb)) = & 
     &                 unk_e_z(:,is_edgez:ie_edgez,js_edgez:je_edgez,  &
     &                           ks_edgez:ke_edgez,point_to(lb))
                  end if
               end if

               if (nvarcorn.gt.0) then
                  unk_n(:,is_unkn:ie_unkn,js_unkn:je_unkn,  &
     &                    ks_unkn:ke_unkn,new_loc(1,lb)) =  & 
     &                 unk_n(:,is_unkn:ie_unkn,js_unkn:je_unkn,  &
     &                         ks_unkn:ke_unkn,point_to(lb))
               end if

               test(point_to(lb)) = -1
             else
               if (nvar.gt.0) then
        unk(:,is_unk:ie_unk,js_unk:je_unk,ks_unk:ke_unk,new_loc(1,lb)) =  & 
     &   unk(:,is_unk:ie_unk,js_unk:je_unk,ks_unk:ke_unk,lb)
               end if
! move facevars
               if (nfacevar.gt.0) then
                  facevarx(:,is_facex:ie_facex,js_facex:je_facex,  &
     &                       ks_facex:ke_facex,new_loc(1,lb)) =  & 
     &                 facevarx(:,is_facex:ie_facex,js_facex:je_facex,  &
     &                            ks_facex:ke_facex,lb)
                  if (ndim >= 2) then
                  facevary(:,is_facey:ie_facey,js_facey:je_facey,  &
     &                       ks_facey:ke_facey,new_loc(1,lb)) =  & 
     &                 facevary(:,is_facey:ie_facey,js_facey:je_facey,  &
     &                            ks_facey:ke_facey,lb)
                  end if
                  if (ndim == 3) then
                  facevarz(:,is_facez:ie_facez,js_facez:je_facez,  &
     &                       ks_facez:ke_facez,new_loc(1,lb)) =  & 
     &                 facevarz(:,is_facez:ie_facez,js_facez:je_facez,  &
     &                            ks_facez:ke_facez,lb)
                  end if
               end if
! move edgevars
               if (nvaredge.gt.0) then
                  unk_e_x(:,is_edgex:ie_edgex,js_edgex:je_edgex,  &
     &                      ks_edgex:ke_edgex,new_loc(1,lb)) = & 
     &                 unk_e_x(:,is_edgex:ie_edgex,js_edgex:je_edgex,  &
     &                           ks_edgex:ke_edgex,lb)
                  if (ndim >= 2) then
                  unk_e_y(:,is_edgey:ie_edgey,js_edgey:je_edgey,  &
     &                      ks_edgey:ke_edgey,new_loc(1,lb)) = & 
     &                 unk_e_y(:,is_edgey:ie_edgey,js_edgey:je_edgey,  &
     &                           ks_edgey:ke_edgey,lb)
                  end if
                  if (ndim == 3) then
                  unk_e_z(:,is_edgez:ie_edgez,js_edgez:je_edgez,  &
     &                      ks_edgez:ke_edgez,new_loc(1,lb)) = & 
     &                 unk_e_z(:,is_edgez:ie_edgez,js_edgez:je_edgez,  &
     &                           ks_edgez:ke_edgez,lb)
                  end if
               end if
               if (nvarcorn.gt.0) then
               unk_n(:,is_unkn:ie_unkn,js_unkn:je_unkn,  &
     &                  ks_unkn:ke_unkn,new_loc(1,lb)) =   &
     &           unk_n(:,is_unkn:ie_unkn,js_unkn:je_unkn,  &
     &                   ks_unkn:ke_unkn,lb)
               end if

               free(lb) = .true.
             end if
             moved(lb) = .true.
            end if
           end if
         end if

         if (lb.le.new_lnblocks) then
            if (free(lb).and..not.sent(lb)) then
               sent(lb) = .true.
               if (.not.newchild(lb)) then
                  if (old_loc(2,lb).ne.mype) then
                     nsend = nsend + 1
                     call MPI_ISEND (free(lb), & 
     &                               1, & 
     &                               MPI_LOGICAL, & 
     &                               old_loc(2,lb), & 
     &                               maxblocks+lb, & 
     &                               amr_mpi_meshComm, & 
     &                               reqs(nsend), & 
     &                               ierr)
                  end if
               end if
            end if
         end if

         if (lb.le.lnblocks_old.and..not.free(lb)) then
            nmoved = nmoved + 1
            point_to(lb) = max(lnblocks_old,new_lnblocks)+nmoved
            if (point_to(lb).gt.maxblocks) then
               do lb2 = max(lnblocks_old,new_lnblocks)+1,maxblocks
                  if (test(lb2).eq.-1) then
                     point_to(lb) = lb2
                     go to 22
                  end if
               end do
            end if
 22         if (point_to(lb).le.maxblocks) then
               test(point_to(lb)) = 1
               if (nvar.gt.0) then
               unk(:,is_unk:ie_unk,js_unk:je_unk,ks_unk:ke_unk,point_to(lb)) =  & 
     &          unk(:,is_unk:ie_unk,js_unk:je_unk,ks_unk:ke_unk,lb)
               end if
! move facevars
               if (nfacevar.gt.0) then
                  facevarx(:,is_facex:ie_facex,js_facex:je_facex,  &
     &                       ks_facex:ke_facex,point_to(lb)) =  & 
     &                 facevarx(:,is_facex:ie_facex,js_facex:je_facex, &
     &                            ks_facex:ke_facex,lb)
                  if (ndim >= 2) then
                  facevary(:,is_facey:ie_facey,js_facey:je_facey,  &
     &                       ks_facey:ke_facey,point_to(lb)) =  & 
     &                 facevary(:,is_facey:ie_facey,js_facey:je_facey,  &
     &                            ks_facey:ke_facey,lb)
                  end if
                  if (ndim == 3) then
                  facevarz(:,is_facez:ie_facez,js_facez:je_facez,  &
     &                       ks_facez:ke_facez,point_to(lb)) =  & 
     &                 facevarz(:,is_facez:ie_facez,js_facez:je_facez,  &
     &                            ks_facez:ke_facez,lb)
                  end if
               end if
! move edgevars
               if (nvaredge.gt.0) then
                  unk_e_x(:,is_edgex:ie_edgex,js_edgex:je_edgex,  &
     &                      ks_edgex:ke_edgex,point_to(lb)) = & 
     &                 unk_e_x(:,is_edgex:ie_edgex,js_edgex:je_edgex,  &
     &                           ks_edgex:ke_edgex,lb)
                  if (ndim >= 2) then
                  unk_e_y(:,is_edgey:ie_edgey,js_edgey:je_edgey,  &
     &                      ks_edgey:ke_edgey,point_to(lb)) = & 
     &                 unk_e_y(:,is_edgey:ie_edgey,js_edgey:je_edgey,  &
     &                           ks_edgey:ke_edgey,lb)
                  end if
                  if (ndim == 3) then
                  unk_e_z(:,is_edgez:ie_edgez,js_edgez:je_edgez,  &
     &                      ks_edgez:ke_edgez,point_to(lb)) = & 
     &                 unk_e_z(:,is_edgez:ie_edgez,js_edgez:je_edgez,  &
     &                           ks_edgez:ke_edgez,lb)
                  end if
               end if
               if (nvarcorn.gt.0) then
               unk_n(:,is_unkn:ie_unkn,js_unkn:je_unkn,  &
     &                 ks_unkn:ke_unkn,point_to(lb)) =   &
     &           unk_n(:,is_unkn:ie_unkn,js_unkn:je_unkn,  &
     &                   ks_unkn:ke_unkn,lb)
               end if
               free(lb) = .TRUE.
            end if
         end if

         return

      else

         return

      end if

      return
      end subroutine send_block_data
