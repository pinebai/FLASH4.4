!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/mpi_amr_get_remote_block_fvar
!! NAME
!!
!!   mpi_amr_get_remote_block_fvar
!!
!! SYNOPSIS
!!
!!   call mpi_amr_get_remote_block_fvar (mype, remote_pe, remote_block, icoord,
!!                                       recvx, recvy, recvz, idest)
!!
!!   call mpi_amr_get_remote_block_fvar (integer, integer, integer, integer,
!!                                       real, real, real, integer)
!!
!! ARGUMENTS
!!
!!   integer, intent(in) :: mype             
!!     The local processor
!!
!!   integer, intent(in) :: remote_pe        
!!     The remote processor.
!!
!!   integer, intent(in) :: remote_block     
!!     The local block id of the block to be copied from
!!     the remote processor.
!!    
!!   integer, intent(in) :: icoord           
!!     Coordinate to fetch, ie facevarx, facveary or facevarz.
!!
!!   real, intent(out) :: recvx            
!!     Output array if icoord = 1.
!!  
!!   real, intent(out) :: recvy
!!     Output array if icoord = 2.
!!
!!   real, intent(out) :: recvz
!!     Output array if icoord = 3.
!!
!!   integer, intent(in) :: idest            
!!     Selects the storage space in the 1blk data structures which is to
!!     be used in this call. If the leaf node is having its
!!     guardcells filled then set this to 1, if its parent
!!     is being filled set it to 2.
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
!!   workspace
!!   mpi_morton
!!   paramesh_interfaces
!!   paramesh_mpi_interfaces
!!
!! CALLS
!!
!!   amr_mpi_find_blk_in_buffer
!!   mpi_set_message_limits
!!
!! RETURNS
!!
!!   Upon return the data from 'remote_block' on 'remote_pe' is placed locally on the
!!   calling processor in the arrays 'recvx', 'recvy', or 'recvz'.
!!
!! DESCRIPTION
!! 
!!  This routine copies guard cell information to face iface in layer
!!  idest of the working block, from the appropriate face of the neighboring 
!!  block, assuming that the neighboring block is on a different processor.
!! 
!! AUTHORS
!!
!! Written by Peter MacNeice (August 2001).
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"


      subroutine mpi_amr_get_remote_block_fvar(mype, & 
     &                     remote_pe,remote_block,icoord, & 
     &                     recvx,recvy,recvz,idest)



      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace
      use mpi_morton
      use paramesh_comm_data

      use paramesh_interfaces, only : amr_mpi_find_blk_in_buffer
      use paramesh_mpi_interfaces, only : mpi_set_message_limits

      implicit none

      include 'mpif.h'


#ifdef TIMINGS
#include "timer.fh"
#endif

!-------------------------
      integer, intent(in) :: mype,remote_pe,remote_block
      integer, intent(in) :: icoord,idest
      real, intent(out)   :: recvx(:,:,:,:)
      real, intent(out)   :: recvy(:,:,:,:)
      real, intent(out)   :: recvz(:,:,:,:)

      integer :: nguard0
      integer :: nguard_work0

      logical :: lfound,lerror
      integer :: ierrorcode,ierr
      integer :: dtype
      integer :: vtype
      integer :: index,index0
      integer :: i, j, k, ia, ib, ja, jb, ka, kb
!-------------------------

      nguard0 = nguard*npgs
      nguard_work0 = nguard_work*npgs

#ifdef TIMINGS
      itimer1 = irtc()
#endif


      if(remote_block.le.lnblocks.and.remote_pe.eq.mype) then


! Copy complete remote block into a buffer block called recv.
        if (no_permanent_guardcells) then
        if(icoord.eq.1) then
          recvx(:,:,:,:) = gt_facevarx(:,:,:,:,remote_block)
        elseif(icoord.eq.2) then
          if (ndim >= 2) then
             recvy(:,:,:,:) = gt_facevary(:,:,:,:,remote_block)
          end if
        elseif(icoord.eq.3) then
          if (ndim == 3) then
             recvz(:,:,:,:) = gt_facevarz(:,:,:,:,remote_block)
          end if
        endif
        else ! no_permanent_guardcells
        if(icoord.eq.1) then
          recvx(:,:,:,:) = facevarx(:,:,:,:,remote_block)
        elseif(icoord.eq.2) then
          if (ndim >= 2) then
             recvy(:,:,:,:) = facevary(:,:,:,:,remote_block)
          end if
        elseif(icoord.eq.3) then
          if (ndim == 3) then
             recvz(:,:,:,:) = facevarz(:,:,:,:,remote_block)
          end if
        endif

        endif ! no_permanent_guardcells

!-----


      else


        call amr_mpi_find_blk_in_buffer(mype,remote_block, & 
     &                        remote_pe,idest,dtype,index0,lfound)

        lerror=.false.
        if(.not.lfound) lerror=.true.
        if( icoord.eq.1 .and. & 
     &      (dtype.ne.13.and.dtype.ne.14.and.dtype.ne.15.and. & 
     &       dtype.ne.13+27.and.dtype.ne.14+27.and.dtype.ne.15+27) & 
     &      ) lerror = .true.
        if( icoord.eq.2 .and. & 
     &      (dtype.ne.11.and.dtype.ne.14.and.dtype.ne.17.and. & 
     &       dtype.ne.11+27.and.dtype.ne.14+27.and.dtype.ne.17+27) & 
     &      ) lerror = .true.
        if( icoord.eq.3 .and. & 
     &      (dtype.ne.5.and.dtype.ne.14.and.dtype.ne.23.and. & 
     &       dtype.ne.5+27.and.dtype.ne.14+27.and.dtype.ne.23+27) & 
     &      ) lerror = .true.

        if(lerror) then
          write(*,*) 'Paramesh error : pe ',mype, & 
     &      ' needed remote blk facevar data ', & 
     &      remote_block,remote_pe,' but could not find it or only ', & 
     &      ' found part of it in the message buffer.', & 
     &      '  Contact developers for help.', & 
     &      ' lfound ',lfound,' dtype ',dtype,' icoord ',icoord
          call mpi_abort(amr_mpi_meshComm,ierrorcode,ierr)
        endif



!-----

! starting index if cell-centered data is also included in recv_buf
        index = index0
        if(l_datapacked(2)) index = & 
     &                      index + nvar*message_size_cc(dtype)

        vtype = 2
        call mpi_set_message_limits( & 
     &               dtype,ia,ib,ja,jb,ka,kb,vtype)

        if(icoord.eq.1) then
          do k = ka,kb
          do j = ja,jb
          do i = ia,ib
            recvx(1:nbndvar,i,j,k) = & 
     &              temprecv_buf(index+1:index+nbndvar)
            index = index+nbndvar
          enddo
          enddo
          enddo
        else
          index  = index + nbndvar*(ib-ia+1)*(jb-ja+1)*(kb-ka+1)
        endif

        if(ndim.ge.2) then
          vtype = 3
          call mpi_set_message_limits( & 
     &               dtype,ia,ib,ja,jb,ka,kb,vtype)

          if(icoord.eq.2) then
            do k = ka,kb
            do j = ja,jb
            do i = ia,ib
              recvy(1:nbndvar,i,j,k) = & 
     &              temprecv_buf(index+1:index+nbndvar)
              index = index+nbndvar
            enddo
            enddo
            enddo
          else
            index  = index + nbndvar*(ib-ia+1)*(jb-ja+1)*(kb-ka+1)
          endif

        endif

        if(ndim.eq.3) then
         vtype = 4
         call mpi_set_message_limits( & 
     &               dtype,ia,ib,ja,jb,ka,kb,vtype)

         if(icoord.eq.3) then
           do k = ka,kb
           do j = ja,jb
           do i = ia,ib
             recvz(1:nbndvar,i,j,k) = & 
     &              temprecv_buf(index+1:index+nbndvar)
             index = index+nbndvar
           enddo
           enddo
           enddo
         endif
        else
         index  = index + nbndvar*(ib-ia+1)*(jb-ja+1)*(kb-ka+1)
        endif

!-----

      endif


#ifdef TIMINGS
      itimer2 = irtc()
      irtc_cprem = itimer2-itimer1+irtc_cprem
#endif

      return
      end subroutine mpi_amr_get_remote_block_fvar
