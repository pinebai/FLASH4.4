!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/amr_edge_diagonal_check
!! NAME
!!
!!   amr_edge_diagonal_check (mype)
!!
!! SYNOPSIS
!!
!!   call amr_edge_diagonal_check (mype)
!!
!!   call amr_edge_diagonal_check (integer)
!!
!! ARGUMENTS
!!
!!   integer, intent(in) :: mype
!!     The calling processors id.
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
!!   paramesh_interfaces
!!   paramesh_mpi_interfaces
!!
!! CALLS
!!
!!   amr_mpi_find_blk_in_buffer
!!   mpi_set_message_limits
!!   mpi_put_edge_buffer_1blk
!!
!! RETURNS
!!
!!   Nothing returned. 
!!
!! DESCRIPTION
!! 
!!  This routine checks to see if the diagonal block between two
!!  leaf-neighbors at the same refinement level as the current block,
!!  is refined. If it is then the edge-based variables along the edge
!!  shared with that diagonal block is given the edge values
!!  form the refined diagonal block, to insure conservation properties.
!! 
!! AUTHORS
!!
!!  Written by Peter MacNeice (October 1997).
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

!#define DEBUG

      subroutine amr_edge_diagonal_check(mype)

      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton

      use paramesh_interfaces, only : amr_mpi_find_blk_in_buffer
      use paramesh_mpi_interfaces, only : mpi_set_message_limits, & 
     &                                    mpi_put_edge_buffer_1blk

      implicit none 

      include 'mpif.h'

      integer, intent(in)  ::  mype

!------------------------------------
! local variables

      integer :: nguard0
      integer :: klo,kup
      integer :: jlo,jup
      integer :: ilo,iup

      integer :: remote_pe,remote_block
      integer :: remote_pe2,remote_block2
      integer :: mark_edge(12,maxblocks)
      integer :: i, ie, iblk, lb, k, j, ierrcode, ierr
      integer :: ia, ib, ja, jb, ka, kb
      integer :: ia0, ib0, ja0, jb0, ka0, kb0
      integer :: dtype, vtype
      integer :: index, index0, n

      logical :: lfound

      real,allocatable :: receive(:,:)

!------------------------------------

      nguard0 = nguard*npgs
      klo=1+nguard0*k3d
      kup=klo*k3d+nzb
      jlo=1+nguard0*k2d
      jup=jlo*k2d+nyb
      ilo=1+nguard0
      iup=ilo+nxb

      allocate(receive(nedges,maxdim+2*nguard0))

#ifdef DEBUG
       call amr_flush(6)
       call mpi_barrier (amr_mpi_meshComm, ierrcode)
       write(*,*) 'starting amr_edge_diagonal_check : pe ',mype
       call amr_flush(6)
       call mpi_barrier (amr_mpi_meshComm, ierrcode)
#endif /* DEBUG */

       if (ndim >= 2) then


! Initialize array marking edges for diagonal patching.
      mark_edge(:,:) = 0

      do i = 1,no_of_diagonal_edges
        ie   = edge_mark(6,1,i)
        iblk = edge_mark(6,2,i)
        mark_edge(ie,iblk) = i
      enddo


! Loop over the blocks on this processor.
      if(lnblocks.gt.0) then
      do lb=1,lnblocks

! Is this a leaf block which has finished its current timestep?
      if((nodetype(lb).eq.1.and..not.var_dt) .or. & 
     &   (nodetype(lb).eq.1.and.ldtcomplete(lb))) then

#ifdef DEBUG
         write(*,*) 'amr_edge_diagonal_check : ', & 
     &              'checking edges on blk ',lb
         write(*,*) 'nbedges ',nbedges
#endif /* DEBUG */


! Any edges on this block which are still marked need a diagonal patch.
! Note that in the data copies below, we can always assume that a
! neighbor block exists, since the edge would not have been marked
! earlier if that was not so.

! Loop over the edges on this block.
       do ie=1,nbedges


       if(mark_edge(ie,lb).ge.1) then

        lfound = .false.
        remote_block  = edge_mark(6,3,mark_edge(ie,lb))
        remote_pe     = edge_mark(6,4,mark_edge(ie,lb))
        remote_block2 = edge_mark(6,3,mark_edge(ie,lb))
        remote_pe2    = edge_mark(6,4,mark_edge(ie,lb))

#ifdef DEBUG
         write(*,*) 'amr_edge_diagonal_check : ', & 
     &              'data source edge for ',lb,  & 
     &              ' is on blk ',remote_block,remote_pe
         write(*,*) 'nbedges ',nbedges
#endif /* DEBUG */

! (remote_block,remote_pe) may be a local block or a remote block.


         if(remote_pe2.ne.mype) then

            lfound = .false.
            do iblk = strt_buffer,last_buffer
               if(remote_block2.eq.laddress(1,iblk).and. & 
     &              remote_pe2 .eq.laddress(2,iblk) ) then
                  remote_block2 = iblk
                  remote_pe2    = mype
                  lfound = .true.
               endif
            enddo

         elseif(remote_pe2.eq.mype) then

            lfound = .true.

         endif


! The edge data on the neighboring faces can be assumed to have been averaged
! correctly from the refined diagonal blocks.

      if (remote_pe == mype .and. remote_block <= lnblocks) then

! Now copy over the edge data from one of the neighbors.
       if(ie.eq.1) then                    ! edge: x low edge, y low edge
         do k=klo,kup-k3d
           bedge_facex_z(:,1,1+nguard0*k2d,k,lb) = &
             bedge_facex_z(:,2,jup,k,remote_block)
         enddo
         bedge_facey_z(:,1+nguard0,1,klo:kup-k3d,lb)= & 
     &                bedge_facex_z(:,1,1+nguard0*k2d,klo:kup-k3d,lb)

       elseif(ie.eq.2) then               ! edge: x low edge, y high edge
         do k=klo,kup-k3d
           bedge_facex_z(:,1,k2d+nguard0*k2d+nyb,k,lb) =  & 
              bedge_facex_z(:,2,jlo,k,remote_block)
         enddo
         bedge_facey_z(:,1+nguard0,2,klo:kup-k3d,lb)= & 
     &     bedge_facex_z(:,1,k2d+nguard0*k2d+nyb,klo:kup-k3d,lb)

       elseif(ie.eq.3) then               ! edge: x high edge, y low edge
         do k=klo,kup-k3d
           bedge_facex_z(:,2,1+nguard0*k2d,k,lb) = &
              bedge_facex_z(:,1,jup,k,remote_block)
         enddo
         bedge_facey_z(:,1+nguard0+nxb,1,klo:kup-k3d,lb)= & 
     &     bedge_facex_z(:,2,1+nguard0*k2d,klo:kup-k3d,lb)


       elseif(ie.eq.4) then               ! edge: x high edge, y high edge
         do k=klo,kup-k3d
           bedge_facex_z(:,2,k2d+nguard0*k2d+nyb,k,lb) =  & 
              bedge_facex_z(:,1,jlo,k,remote_block)
         enddo
         bedge_facey_z(:,1+nguard0+nxb,2,klo:kup-k3d,lb)= & 
     &     bedge_facex_z(:,2,k2d+nguard0*k2d+nyb,klo:kup-k3d,lb)


       elseif(ie.eq.5) then                ! edge: y low edge, z low edge
         do i=ilo,iup-1
           bedge_facey_x(:,i,1,klo,lb) =  &
              bedge_facey_x(:,i,2,kup,remote_block)
         enddo
         bedge_facez_x(:,ilo:iup-1,1+nguard0*k3d,1,lb)= & 
     &                bedge_facey_x(:,ilo:iup-1,1,klo,lb)


       elseif(ie.eq.6) then                ! edge: y high edge, z low edge
         do i=ilo,iup-1
           bedge_facey_x(:,i,2,klo,lb) =  &
              bedge_facey_x(:,i,1,kup,remote_block)
         enddo
         bedge_facez_x(:,ilo:iup-1,k2d+nguard0*k2d+nyb,1,lb)= & 
     &                bedge_facey_x(:,ilo:iup-1,2,klo,lb)


       elseif(ie.eq.7) then                ! edge: y low edge, z high edge
         do i=ilo,iup-1
           bedge_facey_x(:,i,1,kup,lb) =  &
              bedge_facey_x(:,i,2,klo,remote_block)
         enddo
         bedge_facez_x(:,ilo:iup-1,1+nguard0*k2d,2,lb)= & 
     &                bedge_facey_x(:,ilo:iup-1,1,kup,lb)


       elseif(ie.eq.8) then                ! edge: y high edge, z high edge
         do i=ilo,iup-1
           bedge_facey_x(:,i,2,kup,lb) =  &
              bedge_facey_x(:,i,1,klo,remote_block)
         enddo
         bedge_facez_x(:,ilo:iup-1,k2d+nguard0*k2d+nyb,2,lb)= & 
     &                bedge_facey_x(:,ilo:iup-1,2,kup,lb)

       elseif(ie.eq.9) then                ! edge: x low edge, z low edge
         do j=jlo,jup-k2d
           bedge_facex_y(:,1,j,klo,lb) =  &
              bedge_facex_y(:,2,j,kup,remote_block)
         enddo
         bedge_facez_y(:,1+nguard0,jlo:jup-k2d,1,lb)= & 
     &                bedge_facex_y(:,1,jlo:jup-k2d,klo,lb)


       elseif(ie.eq.10) then                ! edge: x low edge, z high edge
         do j=jlo,jup-k2d
           bedge_facex_y(:,1,j,kup,lb) =  &
              bedge_facex_y(:,2,j,klo,remote_block)
         enddo
         bedge_facez_y(:,1+nguard0,jlo:jup-k2d,2,lb)= & 
     &                bedge_facex_y(:,1,jlo:jup-k2d,kup,lb)


       elseif(ie.eq.11) then                ! edge: x high edge, z low edge
         do j=jlo,jup-k2d
           bedge_facex_y(:,2,j,klo,lb) =  &
              bedge_facex_y(:,1,j,kup,remote_block)
         enddo
         bedge_facez_y(:,1+nguard0+nxb,jlo:jup-k2d,1,lb)= & 
     &                bedge_facex_y(:,2,jlo:jup-k2d,klo,lb)


       elseif(ie.eq.12) then                ! edge: x high edge, z high edge
         do j=jlo,jup-k2d
           bedge_facex_y(:,2,j,kup,lb) =  &
              bedge_facex_y(:,1,j,klo,remote_block)
         enddo
         bedge_facez_y(:,1+nguard0+nxb,jlo:jup-k2d,2,lb)= & 
     &                bedge_facex_y(:,2,jlo:jup-k2d,kup,lb)
       endif


      else                      ! if (remote_pe

         call mpi_put_edge_buffer_1blk(lb,remote_block,remote_pe)

         if(ie.eq.1) then       ! edge: x low edge, y low edge
            do k=klo,kup-k3d
               bedge_facex_z(:,1,jlo,k,lb)=  & 
     &              recvarx2e(:,2,jup,k)
            enddo
            bedge_facey_z(:,ilo,1,klo:kup-k3d,lb)= & 
     &           bedge_facex_z(:,1,jlo,klo:kup-k3d,lb)
            
         elseif(ie.eq.2) then   ! edge: x low edge, y high edge
            do k=klo,kup-k3d
               bedge_facex_z(:,1,jup,k,lb)=  & 
     &              recvarx2e(:,2,jlo,k)
            enddo
            bedge_facey_z(:,ilo,2,klo:kup-k3d,lb)= & 
     &           bedge_facex_z(:,1,jup,klo:kup-k3d,lb)
            
         elseif(ie.eq.3) then   ! edge: x high edge, y low edge
            do k=klo,kup-k3d
               bedge_facex_z(:,2,jlo,k,lb)=  & 
     &              recvarx2e(:,1,jup,k)
            enddo
            bedge_facey_z(:,iup,1,klo:kup-k3d,lb)= & 
     &           bedge_facex_z(:,2,jlo,klo:kup-k3d,lb)
            
         elseif(ie.eq.4) then   ! edge: x high edge, y high edge
            do k=klo,kup-k3d
               bedge_facex_z(:,2,jup,k,lb)=  & 
     &              recvarx2e(:,1,jlo,k)
            enddo
            bedge_facey_z(:,iup,2,klo:kup-k3d,lb)= & 
     &           bedge_facex_z(:,2,jup,klo:kup-k3d,lb)
         elseif(ie.eq.5) then   ! edge: y low edge, z low edge
            do i=ilo,iup-1
               bedge_facey_x(:,i,1,klo,lb)= recvary1e(:,i,2,kup)
            enddo
            bedge_facez_x(:,ilo:iup-1,jlo,1,lb)= & 
     &           bedge_facey_x(:,ilo:iup-1,1,klo,lb)
            
         elseif(ie.eq.6) then   ! edge: y high edge, z low edge
            do i=ilo,iup-1
               bedge_facey_x(:,i,2,klo,lb)= recvary1e(:,i,1,kup)
            enddo
            bedge_facez_x(:,ilo:iup-1,jup,1,lb)= & 
     &           bedge_facey_x(:,ilo:iup-1,2,klo,lb)
            
         elseif(ie.eq.7) then   ! edge: y low edge, z high edge
            do i=ilo,iup-1
               bedge_facey_x(:,i,1,kup,lb)= recvary1e(:,i,2,klo)
            enddo
            bedge_facez_x(:,ilo:iup-1,jlo,2,lb)= & 
     &           bedge_facey_x(:,ilo:iup-1,1,kup,lb)
            
         elseif(ie.eq.8) then   ! edge: y high edge, z high edge
            do i=ilo,iup-1
               bedge_facey_x(:,i,2,kup,lb)= recvary1e(:,i,1,klo)
            enddo
            bedge_facez_x(:,ilo:iup-1,jup,2,lb)= & 
     &           bedge_facey_x(:,ilo:iup-1,2,kup,lb)
            
         elseif(ie.eq.9) then   ! edge: x low edge, z low edge
            do j=jlo,jup-k2d
               bedge_facex_y(:,1,j,klo,lb)= recvarx1e(:,2,j,kup)
            enddo
            bedge_facez_y(:,1+nguard0,jlo:jup-k2d,1,lb)= & 
     &           bedge_facex_y(:,1,jlo:jup-k2d,klo,lb)
            
         elseif(ie.eq.10) then  ! edge: x low edge, z high edge
            do j=jlo,jup-k2d
               bedge_facex_y(:,1,j,kup,lb)= recvarx1e(:,2,j,klo)
            enddo
            bedge_facez_y(:,1+nguard0,jlo:jup-k2d,2,lb)= & 
     &           bedge_facex_y(:,1,jlo:jup-k2d,kup,lb)
            
         elseif(ie.eq.11) then  ! edge: x high edge, z low edge
            do j=jlo,jup-k2d
               bedge_facex_y(:,2,j,klo,lb)= recvarx1e(:,1,j,kup)
            enddo
            bedge_facez_y(:,1+nguard0+nxb,jlo:jup-k2d,1,lb)= & 
     &           bedge_facex_y(:,2,jlo:jup-k2d,klo,lb)
            
         elseif(ie.eq.12) then  ! edge: x high edge, z high edge
            do j=jlo,jup-k2d
               bedge_facex_y(:,2,j,kup,lb)= recvarx1e(:,1,j,klo)
            enddo
            bedge_facez_y(:,1+nguard0+nxb,jlo:jup-k2d,2,lb)= & 
     &           bedge_facex_y(:,2,jlo:jup-k2d,kup,lb)
            
         end if

      end if                    ! end if (remote_pe

      endif


      enddo                     ! loop over edges

      endif

      enddo
      endif

      end if ! if (ndim >= 2)

      deallocate(receive)

      return
      end subroutine amr_edge_diagonal_check
