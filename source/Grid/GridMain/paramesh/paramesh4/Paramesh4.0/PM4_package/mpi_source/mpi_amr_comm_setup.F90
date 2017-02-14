!!!#define DEBUG
!!!#define AIX
#if N_DIM > 1
#define DEBUG_DAT
#endif
!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/mpi_amr_comm_setup
!! NAME
!!
!!   mpi_amr_comm_setup
!! 
!! SYNOPSIS
!!
!!   call mpi_amr_comm_setup(mype,nprocs,
!!                           lguard,lprolong,
!!                           lflux,ledge,lrestrict,lfulltree,
!!                           iopt,lcc,lfc,lec,lnc,tag_offset,
!!                           nlayersx,nlayersy,nlayersz,
!!                           flux_dir)
!!   call mpi_amr_comm_setup(mype,nprocs,
!!                           lguard,lprolong,
!!                           lflux,ledge,lrestrict,lfulltree,
!!                           iopt,lcc,lfc,lec,lnc,tag_offset)
!!
!!   call mpi_amr_comm_setup(integer, integer,
!!                           logical, logical,
!!                           logical, logical, logical, logical,
!!                           integer, logical, logical, logical, logical, integer
!!                           optional integer, optional integer, optional, integer,
!!                           optional integer)
!!
!! ARGUMENTS      
!!
!!    integer, intent(in) :: mype,nprocs,iopt
!!      mype   -> The local processor number.
!!      nprocs -> The total number of processors being used.
!!       
!!   logical, intent(in)  :: lguard,lprolong,lflux,ledge,lrestrict,lfulltree
!!      lguard    -> If true, info for guardcell filling will be
!!                    communicated.
!!      lprolong  -> If true, info for prolongation will be
!!                    communicated.
!!      lflux     -> If true, info for flux conservation will be
!!                    communicated.
!!      ledge     -> If true, info for edge averaging will be
!!                    communicated.
!!      lrestrict -> If true, info for restriction will be
!!                    communicated.
!!      lfulltree -> If true, info for restricting data up the entire
!!                    tree will be communicated.
!!
!!   integer, intent(in) :: iopt
!!      iopt   -> Controls which arrays are operated on:
!!                If iopt = 1 then 'unk', 'facevarx(y,z)', 'unk_e_x(y,z)'
!!                and 'unk_n'
!!                If iopt >= 2 then 'work'.
!!
!!   logical, intent(in) :: lcc,lfc,lec,lnc
!!      lcc -> A logical switch controlling whether unk or work data
!!              is filled.
!!      lfc -> A logical switch controlling whether facevar data
!!              is filled.
!!      lec -> A logical switch controlling whether unk_e_? data
!!              is filled.
!!      lnc -> A logical switch controlling whether unk_n data
!!              is filled.
!!
!!   integer, intent(inout) :: tag_offset
!!      tag_offset -> Defines the last tag number used for an mpi message.
!!                    This can be almost anything, but 100 is usually a 
!!                    good choice.
!!
!!   integer, intent(in), optional :: nlayersx,nlayersy,nlayersz
!!     Optional integer arguments specifying how many guardcells to 
!!      exchange in each coordinate direction.  If these are not 
!!      specified then all guardcells are filled.
!!
!!   integer, intent(in), optional :: flux_dir
!!     Optional argument specifying which direction to operate on 
!!       for flux conservation. 
!!       If flux_dir = 1, operate on 'x' direction.
!!       If flux_dir = 2, operate on 'y' direction.
!!       If flux_dir = 3, operate on 'z' direction.
!!     If not specified, all directions are operated on.
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
!!   workspace
!!   tree
!!   timings
!!   mpi_morton
!!   paramesh_mpi_interfaces
!!
!! CALLS
!!
!!   mpi_pack_blocks
!!   mpi_Sbuffer_size
!!   mpi_unpack_blocks
!!   mpi_Rbuffer_size
!!   mpi_pack_edges
!!   mpi_unpack_edges
!!   mpi_pack_fluxes
!!   mpi_unpack_fluxes
!!   mpi_pack_tree_info
!!   mpi_unpack_tree_info
!!   mpi_amr_read_guard_comm
!!   mpi_amr_read_prol_comm
!!   mpi_amr_read_flux_comm
!!   mpi_amr_read_restrict_comm
!!   mpi_set_message_sizes
!!   mpi_xchange_blocks
!!
!! RETURNS
!!
!!   Upon exit, all necessary block data has been communicated to the local
!!   calling processor for the operation selected through the specified
!!   arguments.
!!
!! DESCRIPTION
!!
!!   This routine is the pre-analysis operation which precedes a code block 
!!   in which interprocessor communications may be required.  What is does
!!   is set up and communicate the necessary block information for performing
!!   various communications tasks.  Those tasks are:
!!     1) guardcell filling
!!     2) prolongation
!!     3) restriction
!!     4) flux conservation
!!     5) edge averaging
!!
!! AUTHORS
!!
!!   Peter MacNeice (June 2000) with modifications by Kevin Olson for
!!   directional guardcell filling and flux conservation.
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"


      subroutine mpi_amr_comm_setup(mype,nprocs, & 
     &                              lguard,lprolong, & 
     &                              lflux,ledge,lrestrict,lfulltree, & 
     &                              iopt,lcc,lfc,lec,lnc,tag_offset, & 
     &                              nlayersx,nlayersy,nlayersz, & 
     &                              flux_dir)


! Arguments:
!------------------------------------------------------------------------
      use paramesh_dimensions
      use physicaldata
      use workspace
      use tree
      use timings
      use mpi_morton
      use paramesh_comm_data


      use paramesh_mpi_interfaces, only : & 
     &                                mpi_pack_blocks, & 
     &                                mpi_Sbuffer_size, & 
     &                                mpi_unpack_blocks, & 
     &                                mpi_Rbuffer_size, & 
     &                                mpi_pack_edges, & 
     &                                mpi_unpack_edges, & 
     &                                mpi_pack_fluxes, & 
     &                                mpi_unpack_fluxes, & 
     &                                mpi_pack_tree_info, & 
     &                                mpi_unpack_tree_info, & 
     &                                mpi_amr_read_guard_comm, & 
     &                                mpi_amr_read_prol_comm, & 
     &                                mpi_amr_read_flux_comm, & 
     &                                mpi_amr_read_restrict_comm, & 
     &                                mpi_set_message_sizes, & 
     &                                mpi_xchange_blocks
#ifdef DEBUG_DAT
      use Logfile_interface, ONLY: Logfile_stamp
#endif

      implicit none

      integer, intent(in)    :: mype,nprocs,iopt
      integer, intent(inout) :: tag_offset
      logical, intent(in)    :: lcc,lfc,lec,lnc,lfulltree
      logical, intent(in)    :: lguard,lprolong,lflux,ledge,lrestrict
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz
      integer, intent(in), optional :: flux_dir

      include 'mpif.h'

      real, save, dimension (:), allocatable :: send_buf
      real, save, dimension (:), allocatable :: recv_buf
      integer :: buffer_dim_send, buffer_dim_recv
      integer :: buffer_dim, buf_dim_bytes
      integer :: bufdim
      integer :: max_blks_sent
      integer :: len_surr_blks
      integer :: itemp
      integer :: offset_tree
      integer :: nlayerstx, nlayersty, nlayerstz
      integer :: flux_dirt
      integer :: ierror

#ifdef DEBUG_DAT
      character(len=32), dimension(2,2) :: block_buff
      character(len=32)                 :: int_to_str
      integer, save :: lastWritten_buffer_dim_send = -1
      integer, save :: lastWritten_buffer_dim_recv = -1
#endif

      double precision :: time1
      double precision :: time2

      integer :: ii,jj
      logical, save :: first_call = .true.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      if (timing_mpi) then
         time1 = mpi_wtime()
         time2 = mpi_wtime()
      endif

#ifdef AIX
      buffer_dim = nvar*iu_bnd1*ju_bnd1*ku_bnd1*maxblocks
      buffer_dim_send = buffer_dim
      buffer_dim_recv = buffer_dim
      if (first_call) then
         allocate(send_buf(buffer_dim))
         allocate(temprecv_buf(buffer_dim))
         allocate(recv_buf(buffer_dim))
         first_call = .false.
#ifdef DEBUG_DAT
         if (mype .eq. 0) then  ! Only do this for MASTER_PE
            write (block_buff(1,1), '(a)') 'buffer_dim_send'
            write (int_to_str, '(i11,a1)') buffer_dim_send, ','
            write (block_buff(1,2),'(a)')  trim(adjustl(int_to_str))

            write (block_buff(2,1), '(a)') 'buffer_dim'
            write (int_to_str, '(i11)') buffer_dim
            write (block_buff(2,2), '(a)') trim(adjustl(int_to_str))

            call Logfile_stamp( block_buff, 2, 2, & 
                 &              '[mpi_amr_comm_setup]')
         end if
#endif
      end if
#endif


#ifdef DEBUG
      write(*,*) 'pe ',mype,' entered mpi_amr_comm_setup: ' & 
     &           ,' max_blks_sent ', & 
     &           max_blks_sent,' tag_offset ',tag_offset, & 
     &           ' nprocs ',nprocs, & 
     &           '  gcell_on_cc ', gcell_on_cc,' iopt ',iopt
#endif /* DEBUG */

      if (iopt == 1) then


! install user selection for guardcell variables on reset defaults.
      if(lguard) then
        int_gcell_on_cc = gcell_on_cc
        int_gcell_on_fc = gcell_on_fc
        int_gcell_on_ec = gcell_on_ec
        int_gcell_on_nc = gcell_on_nc
        lguard_in_progress = .true.

        if(nvar.gt.0) then
          ngcell_on_cc = 0
          do ii = 1,nvar
            if(int_gcell_on_cc(ii)) then
              ngcell_on_cc =  ngcell_on_cc + 1
              gcell_on_cc_pointer(ngcell_on_cc) = ii
            endif
          enddo
        endif

        if(nfacevar.gt.0) then
          ngcell_on_fc = 0
          do ii = 1,nfacevar
          do jj = 1,3
            if(int_gcell_on_fc(jj,ii)) then
              ngcell_on_fc(jj) =  ngcell_on_fc(jj) + 1
              gcell_on_fc_pointer(jj,ngcell_on_fc(jj)) = ii
            endif
          enddo
          enddo
        endif

        if(nvaredge.gt.0) then
          ngcell_on_ec = 0
          do ii = 1,nvaredge
          do jj = 1,3
            if(int_gcell_on_ec(jj,ii)) then
              ngcell_on_ec(jj) =  ngcell_on_ec(jj) + 1
              gcell_on_ec_pointer(jj,ngcell_on_ec(jj)) = ii
            endif
          enddo
          enddo
        endif

        if(nvarcorn.gt.0) then
          ngcell_on_nc = 0
          do ii = 1,nvarcorn
            if(int_gcell_on_nc(ii)) then
              ngcell_on_nc =  ngcell_on_nc + 1
              gcell_on_nc_pointer(ngcell_on_nc) = ii
            endif
          enddo
        endif

      else                                    ! lguard

        int_gcell_on_cc(:) = .true.
        int_gcell_on_fc(:,:) = .true.
        int_gcell_on_ec(:,:) = .true.
        int_gcell_on_nc(:) = .true.
        lguard_in_progress = .false.
        ngcell_on_cc = nvar
        do ii=1,nvar
          gcell_on_cc_pointer(ii) = ii
        enddo
        ngcell_on_fc = nfacevar
        do ii=1,nfacevar
          gcell_on_fc_pointer(:,ii) = ii
        enddo
        ngcell_on_ec = nvaredge
        do ii=1,nvaredge
          gcell_on_ec_pointer(:,ii) = ii
        enddo
        ngcell_on_nc = nvarcorn
        do ii=1,nvarcorn
          gcell_on_nc_pointer(ii) = ii
        enddo

      endif                                    ! lguard

      if (present(nlayersx)) then
      nlayerstx = nlayersx
      else
      nlayerstx = nguard
      end if

      if (present(nlayersy)) then
      nlayersty = nlayersy
      else
      nlayersty = nguard
      end if

      if (present(nlayersz)) then
      nlayerstz = nlayersz
      else
      nlayerstz = nguard
      end if

      else

      if (present(nlayersx)) then
      nlayerstx = nlayersx
      else
      nlayerstx = nguard_work
      end if

      if (present(nlayersy)) then
      nlayersty = nlayersy
      else
      nlayersty = nguard_work
      end if

      if (present(nlayersz)) then
      nlayerstz = nlayersz
      else
      nlayerstz = nguard_work
      end if

      endif ! end of if (iopt

      if (lguard) then
         if (nxb/nguard < 2) nlayerstx = min(nlayerstx+1,   nguard)
         if (nyb/nguard < 2) nlayersty = min(nlayersty+k2d, nguard)
         if (nzb/nguard < 2) nlayerstz = min(nlayerstz+k3d, nguard)
      end if
      
      if (present(flux_dir)) then
         flux_dirt = flux_dir
      else
         flux_dirt = 0
      end if

      call mpi_set_message_sizes(iopt, & 
     &     nlayerstx,nlayersty,nlayerstz)

      if (timing_mpi) then
      timer_amr_comm_setup(1) =  timer_amr_comm_setup(1) & 
     &                          + mpi_wtime() - time2
      time2 = mpi_wtime()
      endif

      if(lguard.and.(.not.lrestrict) & 
     &   .or. lfulltree ) then

#ifdef DEBUG
      write(*,*) 'pe ',mype,' calling mpi_amr_read_guard_comm'
#endif /* DEBUG */

         call mpi_amr_read_guard_comm(nprocs)

#ifdef DEBUG
      write(*,*) 'pe ',mype,' exited mpi_amr_read_guard_comm'
#endif /* DEBUG */

      elseif(lprolong) then

         call mpi_amr_read_prol_comm(nprocs)

      elseif((lflux.or.ledge).and.(.not.lrestrict)) then

#ifdef DEBUG
      write(*,*) 'pe ',mype,' calling mpi_amr_read_flux_comm'
#endif /* DEBUG */

        call mpi_amr_read_flux_comm(nprocs)

      elseif(lrestrict) then

#ifdef DEBUG
      write(*,*) 'pe ',mype,' calling mpi_amr_read_restrict_comm'
#endif /* DEBUG */

        call mpi_amr_read_restrict_comm(nprocs)

#ifdef DEBUG
      write(*,*) 'pe ',mype,' exited mpi_amr_read_restrict_comm'
#endif /* DEBUG */

      endif

      if (timing_mpi) then
      timer_amr_comm_setup(2) =  timer_amr_comm_setup(2) & 
     &                          + mpi_wtime() - time2
      time2 = mpi_wtime()
      endif

      itemp = max(sum(commatrix_send), sum(commatrix_recv))
      call MPI_ALLREDUCE (itemp, & 
     &                    max_blks_sent, & 
     &                    1, & 
     &                    MPI_INTEGER, & 
     &                    MPI_MAX, & 
     &                    amr_mpi_meshComm, & 
     &                    ierror)

! If lrestrict is true then we only need tree information.
      if( (.not.lrestrict) .or. (lrestrict.and.lflux) & 
     &                     .or. (lrestrict.and.ledge) & 
     &                     .or. (lrestrict.and.lguard) & 
     &                                            ) then

#ifndef AIX
      if(lguard.or.lprolong) then

        call mpi_Sbuffer_size(mype,nprocs,iopt,lcc,lfc,lec,lnc, & 
     &                        buffer_dim_send,offset_tree, & 
     &                        .true., .false., .false., flux_dir, &
     &                        nlayerstx,nlayersty,nlayerstz)

        call mpi_Rbuffer_size(mype,nprocs,iopt,lcc,lfc,lec,lnc, & 
     &                        buffer_dim_recv, & 
     &                        .true.,.false., .false., flux_dir, &
     &                        nlayerstx,nlayersty,nlayerstz)

      elseif (lflux) then	

        call mpi_Sbuffer_size(mype,nprocs,iopt,lcc,lfc,lec,lnc, & 
     &                        buffer_dim_send,offset_tree, & 
     &                        .false., .true., .false., flux_dir, &
     &                        nlayerstx,nlayersty,nlayerstz)

        call mpi_Rbuffer_size(mype,nprocs,iopt,lcc,lfc,lec,lnc, & 
     &                        buffer_dim_recv, & 
     &                        .false.,.true., .false., flux_dir, &
     &                        nlayerstx,nlayersty,nlayerstz)

      elseif (ledge) then	

        call mpi_Sbuffer_size(mype,nprocs,iopt,lcc,lfc,lec,lnc, & 
     &                        buffer_dim_send,offset_tree, & 
     &                        .false., .false., .true., flux_dir, &
     &                        nlayerstx,nlayersty,nlayerstz)

        call mpi_Rbuffer_size(mype,nprocs,iopt,lcc,lfc,lec,lnc, & 
     &                        buffer_dim_recv, & 
     &                        .false.,.false., .true., flux_dir, &
     &                        nlayerstx,nlayersty,nlayerstz)

      end if	

#ifdef DEBUG_DAT
      if (lguard .and. .not. lrestrict) then
         if (mype .eq. 0) then  ! Only do this for MASTER_PE
            if ((lastWritten_buffer_dim_send .NE. buffer_dim_send) .OR. & 
     &           (lastWritten_buffer_dim_recv   .NE. buffer_dim_recv)) then
               write (block_buff(1,1), '(a)') 'buffer_dim_send'
               write (int_to_str, '(i11,a1)') buffer_dim_send, ','
               write (block_buff(1,2),'(a)')  trim(adjustl(int_to_str))

               write (block_buff(2,1), '(a)') 'buffer_dim_recv'
               write (int_to_str, '(i11)') buffer_dim_recv
               write (block_buff(2,2), '(a)') trim(adjustl(int_to_str))
               
               call Logfile_stamp( block_buff, 2, 2, & 
     &              '[mpi_amr_comm_setup]')
               lastWritten_buffer_dim_send = buffer_dim_send
               lastWritten_buffer_dim_recv = buffer_dim_recv
            end if
         end if
      end if
#endif
#endif


!----
! Set up buffer size info, including space for tree data
#ifdef SAVE_MORTS
      len_surr_blks = 3*3*(1+2*k2d)*(1+2*k3d) & 
     &                + 6*3*(1+2*k2d)*(1+2*k3d)
#else
      len_surr_blks = 3*3*(1+2*k2d)*(1+2*k3d)
#endif
      offset_tree = 32+16+len_surr_blks


#ifdef DEBUG
      write(*,*) 'pe ',mype,' logic switches ', & 
     &                 ' lrestrict,lflux,lguard,ledge,lprolong ', & 
     &                 lrestrict,lflux,lguard,ledge,lprolong
      write(*,*) 'pe ',mype,' setup : buffer_dim ', & 
     &           buffer_dim_send, buffer_dim_recv
#endif /* DEBUG */

#ifndef AIX

      ! changed order of (de/)allocations in next 4 lines - KW 2009-01-30
      if(allocated(send_buf)) deallocate(send_buf)
      if(allocated(temprecv_buf)) deallocate(temprecv_buf)
      allocate( temprecv_buf(buffer_dim_recv))
      allocate( send_buf(buffer_dim_send))

#endif

#ifdef DEBUG
      write(*,*) 'pe ',mype,' setup : send_buf allocated ', & 
     &     ' size ',buffer_dim_send,buffer_dim_recv
#endif /* DEBUG */

      if (timing_mpi) then
      timer_amr_comm_setup(3) =  timer_amr_comm_setup(3) & 
     &                          + mpi_wtime() - time2
      time2 = mpi_wtime()
      endif

      if(lguard.or.lprolong) then

        call mpi_pack_blocks(mype,nprocs,iopt,lcc,lfc,lec,lnc, & 
     &                       buffer_dim_send,send_buf,offset_tree, & 
     &                       nlayerstx,nlayersty,nlayerstz)

      elseif(lflux) then

        call mpi_pack_fluxes(mype,nprocs,buffer_dim_send,send_buf, & 
     &                       offset_tree,flux_dir)

      elseif(ledge) then

        call mpi_pack_edges(mype,nprocs,buffer_dim_send,send_buf, & 
     &                      offset_tree)

      endif

#ifdef DEBUG
      write(*,*) 'pe ',mype,' setup : exited mpi_pack_blocks'
#endif /* DEBUG */

#ifdef DEBUG
      write(*,*) 'pe ',mype,' setup : entering mpi_xchange_blocks' & 
     &       ,' buffer_dims ',buffer_dim_send,buffer_dim_recv
#endif /* DEBUG */
      if (timing_mpi) then
      timer_amr_comm_setup(4) =  timer_amr_comm_setup(4) & 
     &                          + mpi_wtime() - time2
      time2 = mpi_wtime()
      endif


      call mpi_xchange_blocks(mype,nprocs,tag_offset, & 
     &                        buffer_dim_send,send_buf, &
     &                        buffer_dim_recv,temprecv_buf)


      if (timing_mpi) then
      timer_amr_comm_setup(5) =  timer_amr_comm_setup(5) & 
     &                          + mpi_wtime() - time2
      time2 = mpi_wtime()
      endif

#ifdef DEBUG
      write(*,*) 'pe ',mype,' setup : exited mpi_xchange_blocks'
#endif /* DEBUG */

      if(lguard.or.lprolong) then

         call mpi_unpack_blocks(mype,iopt,lcc,lfc,lec,lnc, & 
     &                          buffer_dim_recv,temprecv_buf, & 
     &                          nlayerstx,nlayersty,nlayerstz)

#ifdef DEBUG
      write(*,*) 'pe ',mype,' setup : unpack blocks from buffer'
#endif /* DEBUG */
      elseif(lflux) then
         
         call mpi_unpack_fluxes(mype,buffer_dim_recv,temprecv_buf, & 
     &                          flux_dir)

#ifdef DEBUG
      write(*,*) 'pe ',mype,' setup : exited mpi_unpack_fluxes'
#endif /* DEBUG */

      elseif(ledge) then

         call mpi_unpack_edges(mype,buffer_dim_recv,temprecv_buf)

#ifdef DEBUG
      write(*,*) 'pe ',mype,' setup : exited mpi_unpack_edges'
#endif /* DEBUG */

      endif

      if (timing_mpi) then
      timer_amr_comm_setup(6) =  timer_amr_comm_setup(6) & 
     &                          + mpi_wtime() - time2
      endif

#ifndef AIX
      if(allocated(send_buf)) deallocate(send_buf)
#endif

#ifdef DEBUG
      write(*,*) 'pe ',mype,' setup : deallocated '
#endif /* DEBUG */

!----

      else

!--------------------------------------------------------------
! 
! Establish tree info which is needed about remote blocks

      if (timing_mpi) then
         time2 = mpi_wtime()
      endif


#ifdef DEBUG
      write(*,*)  & 
     &    'pe ',mype,' starting tree setup : max_blks_sent ', & 
     &           max_blks_sent,' offset ',tag_offset
#endif /* DEBUG */

!------
! If we change this buffer size, remember to make the same change in
! mpi_amr_tree_setup
!------
!
! Set up buffer size info for tree messages
#ifdef SAVE_MORTS
      len_surr_blks = 3*3*(1+2*k2d)*(1+2*k3d) & 
     &                + 6*3*(1+2*k2d)*(1+2*k3d)
#else
      len_surr_blks = 3*3*(1+2*k2d)*(1+2*k3d)
#endif

#ifndef AIX
      buffer_dim = (32+16+len_surr_blks)*max_blks_sent

#ifdef DEBUG
      write(*,*) 'pe ',mype,' setup : buffer_dim ', & 
     &           buffer_dim
#endif /* DEBUG */
      
      if(allocated(send_buf)) deallocate(send_buf)
      allocate( send_buf(buffer_dim))
      if(allocated(recv_buf)) deallocate(recv_buf)
      allocate( recv_buf(buffer_dim))
#endif /*AIX*/

#ifdef DEBUG
      write(*,*) 'pe ',mype,' setup : allocated '
#endif /* DEBUG */

      if (timing_mpi) then
      timer_amr_comm_setup(7) =  timer_amr_comm_setup(7) & 
     &                          + mpi_wtime() - time2
      time2 = mpi_wtime()
      endif


! Pack the tree info to be sent

#ifdef DEBUG
      write(*,*) 'pe ',mype,' tree setup : entering mpi_pack_tree_info'
#endif /* DEBUG */
      call mpi_pack_tree_info(mype,nprocs,buf_dim_bytes, & 
     &                     buffer_dim,send_buf)

#ifdef DEBUG
      write(*,*) 'pe ',mype,' tree setup : exited mpi_pack_tree_info'
#endif /* DEBUG */

#ifdef DEBUG
      write(*,*)  & 
     &'pe ',mype,' tree setup : entering mpi_xchange_blocks', & 
     & ' args : mype,nprocs,tag_offset,buffer_dim ', & 
     &  mype,nprocs,tag_offset,buffer_dim
#endif /* DEBUG */

      if (timing_mpi) then
      timer_amr_comm_setup(8) =  timer_amr_comm_setup(8) & 
     &                          + mpi_wtime() - time2
      time2 = mpi_wtime()
      endif

!
! Exchange the tree info to be sent

      call mpi_xchange_blocks(mype,nprocs,tag_offset, & 
     &                        buffer_dim, send_buf, &
     &                        buffer_dim, recv_buf)

      if (timing_mpi) then
      timer_amr_comm_setup(5) =  timer_amr_comm_setup(5) & 
     &                          + mpi_wtime() - time2
      time2 = mpi_wtime()
      endif

#ifdef DEBUG
      write(*,*)  & 
     &'pe ',mype,' tree setup : exited mpi_xchange_blocks'
#endif /* DEBUG */

! Unpack the tree info which has been received
#ifdef DEBUG
      write(*,*)  & 
     & 'pe ',mype,' comm setup : entering mpi_unpack_treeinfo' & 
     &  ,' lguard ',lguard,' lprolong ',lprolong,' lflux ', & 
     &     lflux,' ledge ',' lrestrict ',lrestrict
#endif /* DEBUG */

      call mpi_unpack_tree_info(mype,nprocs,buf_dim_bytes, & 
     &                       buffer_dim,recv_buf)

      if (timing_mpi) then
      timer_amr_comm_setup(9) =  timer_amr_comm_setup(9) & 
     &                          + mpi_wtime() - time2
      time2 = mpi_wtime()
      endif

#ifdef DEBUG
      write(*,*)  & 
     & 'pe ',mype,' tree setup : exited mpi_unpack_treeinfo'
#endif /* DEBUG */

#ifndef AIX
      if(allocated(send_buf)) deallocate(send_buf)
      if(allocated(recv_buf)) deallocate(recv_buf)
#endif

#ifdef DEBUG
      write(*,*) 'pe ',mype,' tree setup : deallocated '
      write(*,*) 'pe ',mype,' exiting mpi_amr_comm_setup ', & 
     &           ' l_datapacked ',l_datapacked
#endif /* DEBUG */

      end if ! end of lrestrict if test


!---------------------------------------------------
! Store laddress pe limits to help optimize searching

      if(last_buffer.gt.strt_buffer.and.nprocs.gt.1) then

         ladd_strt(1:nprocs-1) = last_buffer
         ladd_strt(0         ) = strt_buffer
         do jj= last_buffer,strt_buffer,-1
           ladd_strt(laddress(2,jj)) = jj
         enddo
         do jj= nprocs-2,0,-1
           ladd_strt(jj) = minval(ladd_strt(jj:jj+1))
         enddo

         ladd_end(0:nprocs-2) = strt_buffer
         ladd_end(nprocs-1  ) = last_buffer
         do jj= strt_buffer,last_buffer
           ladd_end(laddress(2,jj)) = jj
         enddo
         do jj= 1,nprocs-1
           ladd_end(jj) = maxval(ladd_end(jj-1:jj))
         enddo

      else
        ladd_strt = strt_buffer
        ladd_end  = strt_buffer
      endif

!---------------------------------------------------


      if (timing_mpi) then
        timer_amr_comm_setup(0) =  timer_amr_comm_setup(0) & 
     &                          + mpi_wtime() - time1
      endif

      return
      end subroutine mpi_amr_comm_setup


