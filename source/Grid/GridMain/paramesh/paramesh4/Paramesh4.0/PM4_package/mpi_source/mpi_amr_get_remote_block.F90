!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/mpi_amr_get_remote_block
!! NAME
!!
!!   mpi_amr_get_remote_block
!!
!! SYNOPSIS
!!
!!   call mpi_amr_get_remote_block (mype, remote_pe, remote_block, idest,
!!                                   iopt, lcc, lfc, lec, lnc,
!!                                   nlayersx,nlayersy,nlayersz )
!!
!!   call mpi_amr_get_remote_block (integer, integer, integer, integer,
!!                                   integer, logical, logical, logical, logical,
!!                                   optional integer, optional integer, 
!!                                   optional integer)
!!
!! ARGUMENTS
!!
!!   integer :: mype             
!!     The local processor
!!
!!   integer :: remote_pe        
!!     The remote processor.
!!
!!   integer :: remote_block     
!!     The local block id of the block to be copied from
!!     the remote processor.
!!    
!!   integer :: idest            
!!     Selects the storage space in the 1blk data structures which is to
!!     be used in this call. If the leaf node is having its
!!     guardcells filled then set this to 1, if its parent
!!     is being filled set it to 2.
!!
!!   integer :: iopt             
!!     A switch to control which data source is to be used
!!     iopt=1 will use 'unk', 'facevarx,y,z', 'unk_e_x,y,z', and 'unk_n'
!!     iopt>=2 will use 'work'
!!
!!   logical :: lcc              
!!     A logical switch which requests cell centered data.
!!
!!   logical :: lfc              
!!     A logical switch which requests cell-face centered data.
!!
!!   logical :: lec              
!!     A logical switch which requests cell-edge centered data.
!!
!!   logical :: lnc              
!!     A logical switch which requests cell-corner data.
!!
!!   optional, interger :: nlayersx, nlayersy, nlayersz
!!     Optional argumentes which indicate the number of guardcells to use
!!     in the x, y, and z directions.
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
!!   Nothing returned. 
!!
!! DESCRIPTION
!! 
!!   This routine copies guard cell information to face iface in layer
!!   idest of the working block, from the appropriate face of the neighboring 
!!   block, assuming that the neighboring block is on a different processor.
!! 
!! AUTHORS
!!
!!  Written by Peter MacNeice (July 1997).  Modified by Kevin Olson for
!!  directional guardcell filling.
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"


      subroutine mpi_amr_get_remote_block(mype,remote_pe,remote_block, & 
     &    idest,iopt,lcc,lfc,lec,lnc, & 
     &    nlayersx,nlayersy,nlayersz)


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
      integer, intent(in) :: idest,iopt
      logical, intent(in) :: lcc,lfc,lec,lnc
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz

      integer :: iopt0


      integer :: nguard0 
      integer :: nguard_work0
      logical :: lfound
      integer :: ierrorcode,ierr
      integer :: dtype
      integer :: vtype
      integer :: index,index0
      integer :: ia, ib, ja, jb, ka, kb, i, j, k
      integer :: ivar, ivar_next


!-------------------------

      nguard0 = nguard*npgs
      nguard_work0 = nguard_work*npgs

#ifdef TIMINGS
      itimer1 = irtc()
#endif


      if(remote_block.le.lnblocks.and.remote_pe.eq.mype) then

!-------------------------
      if(iopt.eq.1) then
!-------------------------

      if(lcc) then

        do ivar=1,nvar
          if(int_gcell_on_cc(ivar)) then
! Copy complete remote block into a buffer block called recv.
          if (no_permanent_guardcells) then
          unk1(ivar, & 
     &         1+nguard:nxb+nguard,1+nguard*k2d:nyb+nguard*k2d, & 
     &         1+nguard*k3d:nzb+nguard*k3d,idest) = & 
     &      gt_unk(ivar, & 
     &         1+nguard0:nxb+nguard0, & 
     &         1+nguard0*k2d:nyb+nguard0*k2d, & 
     &         1+nguard0*k3d:nzb+nguard0*k3d,remote_block)
        else ! no_permanent_guardcells
        unk1(ivar,nguard+1:nguard+nxb, & 
     &              nguard*k2d+1:nguard*k2d+nyb, & 
     &              nguard*k3d+1:nguard*k3d+nzb,idest) = & 
     &   unk(ivar,nguard0+1:nguard0+nxb, & 
     &             nguard0*k2d+1:nguard0*k2d+nyb, & 
     &             nguard0*k3d+1:nguard0*k3d+nzb,remote_block)
         endif ! no_permanent_guardcells
         endif
       enddo                             ! ivar do loop

      endif                              ! lcc if test


!-----

      if(lfc) then

      if(iopt.eq.1) then
      do ivar=1,nbndvar
      if (no_permanent_guardcells) then
      if(int_gcell_on_fc(1,ivar)) then
      facevarx1(ivar, & 
     &          1+nguard:nxb+nguard+1,1+nguard*k2d:nyb+nguard*k2d, & 
     &          1+nguard*k3d:nzb+nguard*k3d,idest) = & 
     & gt_facevarx(ivar, & 
     &             1+nguard0:nxb+nguard0+1, & 
     &             1+nguard0*k2d:nyb+nguard0*k2d, & 
     &             1+nguard0*k3d:nzb+nguard0*k3d,remote_block)
      endif
      if(ndim.ge.2) then
      if(int_gcell_on_fc(2,ivar)) then
      facevary1(ivar, & 
     &          1+nguard:nxb+nguard,1+nguard*k2d:nyb+nguard*k2d+k2d, & 
     &          1+nguard*k3d:nzb+nguard*k3d,idest) = & 
     & gt_facevary(ivar, & 
     &             1+nguard0:nxb+nguard0, & 
     &             1+nguard0*k2d:nyb+nguard0*k2d+k2d, & 
     &             1+nguard0*k3d:nzb+nguard0*k3d,remote_block)
      endif
      endif
      if(ndim.eq.3) then
      if(int_gcell_on_fc(3,ivar)) then
      facevarz1(ivar, & 
     &          1+nguard:nxb+nguard,1+nguard*k2d:nyb+nguard*k2d, & 
     &          1+nguard*k3d:nzb+nguard*k3d+k3d,idest) = & 
     & gt_facevarz(ivar, & 
     &             1+nguard0:nxb+nguard0, & 
     &             1+nguard0*k2d:nyb+nguard0*k2d, & 
     &             1+nguard0*k3d:nzb+nguard0*k3d+k3d,remote_block)
      endif
      end if

      else ! no_permanent_guardcells

      if(int_gcell_on_fc(1,ivar)) then
      facevarx1(ivar, & 
     &          1+nguard:nxb+nguard+1, & 
     &          1+nguard*k2d:nyb+nguard*k2d, & 
     &          1+nguard*k3d:nzb+nguard*k3d,idest) = & 
     &    facevarx(ivar, & 
     &             1+nguard0:nxb+nguard0+1, & 
     &             1+nguard0*k2d:nyb+nguard0*k2d, & 
     &             1+nguard0*k3d:nzb+nguard0*k3d,remote_block)
      endif
      if(ndim.ge.2) then
      if(int_gcell_on_fc(2,ivar)) then
      facevary1(ivar, & 
     &          1+nguard:nxb+nguard,1+nguard*k2d:nyb+nguard*k2d+k2d, & 
     &          1+nguard*k3d:nzb+nguard*k3d,idest) = & 
     &    facevary(ivar, & 
     &             1+nguard0:nxb+nguard0, & 
     &             1+nguard0*k2d:nyb+nguard0*k2d+k2d, & 
     &             1+nguard0*k3d:nzb+nguard0*k3d,remote_block)
      endif
      end if
      if(ndim.eq.3) then
      if(int_gcell_on_fc(3,ivar)) then
      facevarz1(ivar, & 
     &          1+nguard:nxb+nguard,1+nguard*k2d:nyb+nguard*k2d, & 
     &          1+nguard*k3d:nzb+nguard*k3d+k3d,idest) = & 
     &    facevarz(ivar, & 
     &             1+nguard0:nxb+nguard0, & 
     &             1+nguard0*k2d:nyb+nguard0*k2d, & 
     &             1+nguard0*k3d:nzb+nguard0*k3d+k3d,remote_block)
      endif
      end if

      endif ! no_permanent_guardcells

      enddo

      endif

      endif                              ! lfc if test

!-----

!-----

      if (ndim > 1) then
      if(lec) then
      if(iopt.eq.1) then
! Copy complete remote block into a buffer block called recv.
      do ivar=1,nbndvare
      if (no_permanent_guardcells) then
      if(int_gcell_on_ec(1,ivar)) then
      unk_e_x1(ivar, & 
     &         1+nguard:nxb+nguard, & 
     &         1+nguard*k2d:nyb+nguard*k2d+k2d, & 
     &         1+nguard*k3d:nzb+nguard*k3d+k3d,idest) = & 
     & gt_unk_e_x(ivar, & 
     &            1+nguard0:nxb+nguard0, & 
     &            1+nguard0*k2d:nyb+nguard0*k2d+k2d, & 
     &            1+nguard0*k3d:nzb+nguard0*k3d+k3d,remote_block)
      endif
      if(int_gcell_on_ec(2,ivar)) then
      unk_e_y1(ivar, & 
     &         1+nguard:nxb+nguard+1, & 
     &         1+nguard*k2d:nyb+nguard*k2d, & 
     &         1+nguard*k3d:nzb+nguard*k3d+k3d,idest) = & 
     & gt_unk_e_y(ivar, & 
     &            1+nguard0:nxb+nguard0+1, & 
     &            1+nguard0*k2d:nyb+nguard0*k2d, & 
     &            1+nguard0*k3d:nzb+nguard0*k3d+k3d,remote_block)
      endif
      if (ndim == 3) then
      if(int_gcell_on_ec(3,ivar)) then
      unk_e_z1(ivar, & 
     &         1+nguard:nxb+nguard+1, & 
     &         1+nguard*k2d:nyb+nguard*k2d+k2d, & 
     &         1+nguard*k3d:nzb+nguard*k3d,idest) = & 
     & gt_unk_e_z(ivar, & 
     &            1+nguard0:nxb+nguard0+1, & 
     &            1+nguard0*k2d:nyb+nguard0*k2d+k2d, & 
     &            1+nguard0*k3d:nzb+nguard0*k3d,remote_block)
      endif
      end if ! if (ndim == 3)

      else ! no_permanent_guardcells

      if(int_gcell_on_ec(1,ivar)) then
      unk_e_x1(ivar, & 
     &         1+nguard:nxb+nguard, & 
     &         1+nguard*k2d:nyb+nguard*k2d+k2d, & 
     &         1+nguard*k3d:nzb+nguard*k3d+k3d,idest) = & 
     &    unk_e_x(ivar, & 
     &            1+nguard0:nxb+nguard0, & 
     &            1+nguard0*k2d:nyb+nguard0*k2d+k2d, & 
     &            1+nguard0*k3d:nzb+nguard0*k3d+k3d,remote_block)
      endif
      if(int_gcell_on_ec(2,ivar)) then
      unk_e_y1(ivar, & 
     &         1+nguard:nxb+nguard+1, & 
     &         1+nguard*k2d:nyb+nguard*k2d, & 
     &         1+nguard*k3d:nzb+nguard*k3d+k3d,idest) = & 
     &    unk_e_y(ivar, & 
     &            1+nguard0:nxb+nguard0+1, & 
     &            1+nguard0*k2d:nyb+nguard0*k2d, & 
     &            1+nguard0*k3d:nzb+nguard0*k3d+k3d,remote_block)
      endif
      if (ndim == 3) then
      if(int_gcell_on_ec(3,ivar)) then
      unk_e_z1(ivar, & 
     &         1+nguard:nxb+nguard+1, & 
     &         1+nguard*k2d:nyb+nguard*k2d+k2d, & 
     &         1+nguard*k3d:nzb+nguard*k3d,idest) = & 
     &    unk_e_z(ivar, & 
     &            1+nguard0:nxb+nguard0+1, & 
     &            1+nguard0*k2d:nyb+nguard0*k2d+k2d, & 
     &            1+nguard0*k3d:nzb+nguard0*k3d,remote_block)
      endif
      end if ! if (ndim == 3)

      endif  ! no_permanent_guardcells

      enddo
      endif
      endif                              ! lec if test
      end if ! if (ndim > 1
!-----

!-----

      if(lnc) then
      if(iopt.eq.1) then
! Copy complete remote block into a buffer block called recv.
      do ivar=1,nvarcorn
      if(int_gcell_on_nc(ivar)) then
      if (no_permanent_guardcells) then
      unk_n1(ivar, & 
     &       1+nguard:nxb+nguard+1, & 
     &       1+nguard*k2d:nyb+nguard*k2d+k2d, & 
     &       1+nguard*k3d:nzb+nguard*k3d+k3d,idest) = & 
     & gt_unk_n(ivar, & 
     &          1+nguard0:nxb+nguard0+1, & 
     &          1+nguard0*k2d:nyb+nguard0*k2d+k2d, & 
     &          1+nguard0*k3d:nzb+nguard0*k3d+k3d,remote_block)

      else ! no_permanent_guardcells

      unk_n1(ivar, & 
     &       1+nguard:nxb+nguard+1, & 
     &       1+nguard*k2d:nyb+nguard*k2d+k2d, & 
     &       1+nguard*k3d:nzb+nguard*k3d+k3d,idest) = & 
     &    unk_n(ivar, & 
     &          1+nguard0:nxb+nguard0+1, & 
     &          1+nguard0*k2d:nyb+nguard0*k2d+k2d, & 
     &          1+nguard0*k3d:nzb+nguard0*k3d+k3d,remote_block)
      endif ! no_permanent_guardcells
      endif
      enddo
      endif
      endif                              ! lnc if test

!-----


!-------------------------
      elseif(iopt.ge.2) then
!-------------------------
      if(lcc) then

        iopt0 = iopt-1
! Copy complete remote block into a buffer block called recvw.
        work1(1+nguard_work:nxb+nguard_work, & 
     &          1+nguard_work*k2d:nyb+nguard_work*k2d, & 
     &          1+nguard_work*k3d:nzb+nguard_work*k3d,idest) =  & 
     &    work( & 
     &         1+nguard_work0:nxb+nguard_work0, & 
     &         1+nguard_work0*k2d:nyb+nguard_work0*k2d, & 
     &         1+nguard_work0*k3d:nzb+nguard_work0*k3d, & 
     &         remote_block,iopt0)
      endif                              ! lcc if test
!-------------------------
      endif
!-------------------------

      else

      call amr_mpi_find_blk_in_buffer(mype,remote_block, & 
     &                        remote_pe,idest,dtype,index0,lfound)


      if((.not.lfound).or.(dtype.ne.14.and.dtype.ne.14+27)) then
          write(*,*) 'Paramesh error : pe ',mype,' needed full blk ', & 
     &      remote_block,remote_pe,' but could not find it or only ', & 
     &      ' found part of it in the message buffer.', & 
     &      '  Contact PARAMESH developers for help.'
          call mpi_abort(amr_mpi_meshComm,ierrorcode,ierr)
      endif

      if(lcc) then

        if(iopt.eq.1) then
          vtype = 1
          index = index0
          call mpi_set_message_limits( & 
     &                 dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                 nlayersx,nlayersy,nlayersz)

          if (no_permanent_guardcells) then
             ia = ia + nguard
             ib = ib + nguard
             ja = ja + nguard*k2d
             jb = jb + nguard*k2d
             ka = ka + nguard*k3d
             kb = kb + nguard*k3d
          end if

          do k = ka,kb
          do j = ja,jb
          do i = ia,ib
            do ivar=1,ngcell_on_cc
              ivar_next = gcell_on_cc_pointer(ivar)
              unk1(ivar_next,i,j,k,idest) =  & 
     &            temprecv_buf(index+ivar)
            enddo
            index = index+ngcell_on_cc
          enddo
          enddo
          enddo

        elseif(iopt.gt.1) then
          vtype = 0
          index = index0
          call mpi_set_message_limits( & 
     &                 dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                 nlayersx,nlayersy,nlayersz)

          if (no_permanent_guardcells) then
             ia = ia + nguard_work
             ib = ib + nguard_work
             ja = ja + nguard_work*k2d
             jb = jb + nguard_work*k2d
             ka = ka + nguard_work*k3d
             kb = kb + nguard_work*k3d
          end if

          do k = ka,kb
          do j = ja,jb
          do i = ia,ib
            work1(i,j,k,idest) =  temprecv_buf(index+1)
            index = index+1
          enddo
          enddo
          enddo

        endif
      
      endif                              ! lcc if test

!-----

      if(lfc) then

        if(iopt.eq.1) then

! starting index if cell-centered data is also included in recv_buf
        index = index0
        if(l_datapacked(2)) index = & 
     &                      index + ngcell_on_cc*message_size_cc(dtype)

        vtype = 2
        call mpi_set_message_limits( & 
     &               dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &               nlayersx,nlayersy,nlayersz)


        if (no_permanent_guardcells) then
          ia = ia + nguard
          ib = ib + nguard
          ja = ja + nguard*k2d
          jb = jb + nguard*k2d
          ka = ka + nguard*k3d
          kb = kb + nguard*k3d
       end if

        do k = ka,kb
        do j = ja,jb
        do i = ia,ib
            do ivar=1,ngcell_on_fc(1)
              ivar_next = gcell_on_fc_pointer(1,ivar)
              facevarx1(ivar_next,i,j,k,idest) =  & 
     &            temprecv_buf(index+ivar)
            enddo
            index = index+ngcell_on_fc(1)
        enddo
        enddo
        enddo

        if(ndim.ge.2) then
         vtype = 3
         call mpi_set_message_limits( & 
     &              dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &              nlayersx,nlayersy,nlayersz)

         if (no_permanent_guardcells) then
          ia = ia + nguard
          ib = ib + nguard
          ja = ja + nguard*k2d
          jb = jb + nguard*k2d
          ka = ka + nguard*k3d
          kb = kb + nguard*k3d
       end if

         do k = ka,kb
         do j = ja,jb
         do i = ia,ib
            do ivar=1,ngcell_on_fc(2)
              ivar_next = gcell_on_fc_pointer(2,ivar)
              facevary1(ivar_next,i,j,k,idest) =  & 
     &            temprecv_buf(index+ivar)
            enddo
            index = index+ngcell_on_fc(2)
         enddo
         enddo
         enddo
        endif

        if(ndim.eq.3) then
         vtype = 4
         call mpi_set_message_limits( & 
     &               dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &               nlayersx,nlayersy,nlayersz)

         if (no_permanent_guardcells) then
          ia = ia + nguard
          ib = ib + nguard
          ja = ja + nguard*k2d
          jb = jb + nguard*k2d
          ka = ka + nguard*k3d
          kb = kb + nguard*k3d
       end if

         do k = ka,kb
         do j = ja,jb
         do i = ia,ib
            do ivar=1,ngcell_on_fc(3)
              ivar_next = gcell_on_fc_pointer(3,ivar)
              facevarz1(ivar_next,i,j,k,idest) =  & 
     &            temprecv_buf(index+ivar)
            enddo
            index = index+ngcell_on_fc(3)
         enddo
         enddo
         enddo
        endif

        endif

      endif                              ! lfc if test

!-----

!-----
      if (ndim > 1) then
      if(lec) then

        if(iopt.eq.1) then

! starting index if cell-centered data is also included in recv_buf
        index = index0
        If (l_datapacked(2)) index =                                   & 
                             index + ngcell_on_cc*message_size_cc(dtype)
        If(l_datapacked(3)) index =                                    &
                            index + ngcell_on_fc(1) *                  &
                                    message_size_fcx(dtype)            &
                                  + ngcell_on_fc(2) *                  &
                                    message_size_fcy(dtype)            &
                                  + ngcell_on_fc(3) *                  &
                                    message_size_fcz(dtype)   
        vtype = 5
        call mpi_set_message_limits( & 
     &               dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &               nlayersx,nlayersy,nlayersz)

        if (no_permanent_guardcells) then
           ia = ia + nguard
           ib = ib + nguard
           ja = ja + nguard*k2d
           jb = jb + nguard*k2d
           ka = ka + nguard*k3d
           kb = kb + nguard*k3d
        end if

        do k = ka,kb
        do j = ja,jb
        do i = ia,ib
            do ivar=1,ngcell_on_ec(1)
              ivar_next = gcell_on_ec_pointer(1,ivar)
              unk_e_x1(ivar_next,i,j,k,idest) =  & 
     &            temprecv_buf(index+ivar)
            enddo
            index = index+ngcell_on_ec(1)
        enddo
        enddo
        enddo

       if(ndim.ge.2) then
        vtype = 6
        call mpi_set_message_limits( & 
     &               dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &               nlayersx,nlayersy,nlayersz)

        if (no_permanent_guardcells) then
           ia = ia + nguard
           ib = ib + nguard
           ja = ja + nguard*k2d
           jb = jb + nguard*k2d
           ka = ka + nguard*k3d
          kb = kb + nguard*k3d
        endif

        do k = ka,kb
        do j = ja,jb
        do i = ia,ib
            do ivar=1,ngcell_on_ec(2)
              ivar_next = gcell_on_ec_pointer(2,ivar)
              unk_e_y1(ivar_next,i,j,k,idest) =  & 
     &            temprecv_buf(index+ivar)
            enddo
            index = index+ngcell_on_ec(2)
        enddo
        enddo
        enddo
       endif

       if (ndim == 3) then
        vtype =7
        call mpi_set_message_limits( & 
     &               dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &               nlayersx,nlayersy,nlayersz)

        if (no_permanent_guardcells) then
           ia = ia + nguard
           ib = ib + nguard
           ja = ja + nguard*k2d
           jb = jb + nguard*k2d
           ka = ka + nguard*k3d
           kb = kb + nguard*k3d
        endif

        do k = ka,kb
        do j = ja,jb
        do i = ia,ib
            do ivar=1,ngcell_on_ec(3)
              ivar_next = gcell_on_ec_pointer(3,ivar)
              unk_e_z1(ivar_next,i,j,k,idest) =  & 
     &            temprecv_buf(index+ivar)
            enddo
            index = index+ngcell_on_ec(3)
        enddo
        enddo
        enddo
       endif ! if (ndim == 3)

       endif

      endif                              ! lec if test
      end if ! if (ndim > 1) then

!-----

!-----

      if(lnc) then

        if(iopt.eq.1) then

! starting index if cell-centered data is also included in recv_buf
        index = index0
        If (l_datapacked(2)) index =                                   & 
                            index + ngcell_on_cc*message_size_cc(dtype)
        If(l_datapacked(3)) index =                                    &
                            index + ngcell_on_fc(1) *                  &
                                    message_size_fcx(dtype)            &
                                  + ngcell_on_fc(2) *                  &
                                    message_size_fcy(dtype)            &
                                  + ngcell_on_fc(3) *                  &
                                    message_size_fcz(dtype)
           If (l_datapacked(4)) index =                                & 
                            index + maxval(ngcell_on_ec(1:ndim))       & 
                                     *message_size_ec(dtype)

        vtype = 8
        call mpi_set_message_limits( & 
     &               dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &               nlayersx,nlayersy,nlayersz)


        if (no_permanent_guardcells) then
           ia = ia + nguard
           ib = ib + nguard
           ja = ja + nguard*k2d
           jb = jb + nguard*k2d
           ka = ka + nguard*k3d
           kb = kb + nguard*k3d
        endif

        do k = ka,kb
        do j = ja,jb
        do i = ia,ib
            do ivar=1,ngcell_on_nc
              ivar_next = gcell_on_nc_pointer(ivar)
              unk_n1(ivar_next,i,j,k,idest) =  & 
     &            temprecv_buf(index+ivar)
            enddo
            index = index+ngcell_on_nc
        enddo
        enddo
        enddo

        endif

      endif                              ! lnc if test

!-----
      endif

#ifdef TIMINGS
      itimer2 = irtc()
      irtc_cprem = itimer2-itimer1+irtc_cprem
#endif

      return
      end subroutine mpi_amr_get_remote_block
