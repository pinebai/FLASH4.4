!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"


      subroutine amr_perm_to_1blk( lcc,lfc,lec,lnc,lb,pe,iopt,idest)



!------------------------------------------------------------------------
!
! This routine copies data to the 1-block working arrays with guardcells
! from the permanent data arrays, which may or may not have permanent
! guardcells, depending on whether NO_PERMANENT_GUARDCELLS is defined 
! in physicaldata.fh.
!
!
! Arguments :
!      lcc          logical       copies cell centered data if true
!      lfc          logical       copies cell face-centered data if true
!      lec          logical       copies cell edge-centered data if true
!      lnc          logical       copies cell corner data if true
!      lb           integer       block from which data is to be copied
!      pe           integer       processor from which data is to be copied
!      iopt         integer       data structure to be copied
!      idest        integer       sets value for last dimension index
!                                  in the 1-blk data arrays
!
!
! Written :     Peter MacNeice          February 1999
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace
      use mpi_morton
      use paramesh_comm_data
      use paramesh_mpi_interfaces, only : mpi_set_message_limits
      use paramesh_interfaces, only : amr_mpi_find_blk_in_buffer

      implicit none

      include 'mpif.h'

      integer, intent(in) ::  lb,pe,iopt,idest
      logical, intent(in) ::  lcc,lfc,lec,lnc


!------------------------------------
! local variables

      integer :: iopt0
      integer :: nguard0,nguard_work0
      integer :: ia, ib, ja, jb, ka, kb
      integer :: i, j, k, ivar, ivar_next
      integer :: vtype,dtype,rem_blk,rem_pe,mype,ierr
      integer :: index,index0
      logical :: lfound

!------------------------------------

      nguard0 = nguard*npgs
      nguard_work0 = nguard_work*npgs

         if(lb.gt.lnblocks) then
           Call MPI_COMM_RANK(amr_mpi_meshComm, mype, ierr)
           rem_blk = lb
           rem_pe  = mype
           call amr_mpi_find_blk_in_buffer(mype,rem_blk, & 
     &                       rem_pe,idest,dtype,index0,lfound)

           if(.not.lfound) then 
                     write(*,*)  & 
     &              'perm to 1blk reporting blk not found', & 
     &              ' mype=',mype,' looked for ',lb, & 
     &              ' where lnblocks=',lnblocks, & 
     &         ' strt_buffer=',strt_buffer,' last_buffer=',last_buffer, & 
     &         ' laddress ',laddress(:,strt_buffer:last_buffer)
             call amr_abort()    ! remove this abort after testing
           endif

         endif

!
! Put block lb's data into the data_1blk.fh datastructures, with the
! appropriate guardcell padding.
          if(iopt.eq.1) then


          if(lcc) then

          if(lb.le.lnblocks) then

            do ivar=1,nvar
              if(int_gcell_on_cc(ivar)) then
            unk1(ivar,1+nguard:nxb+nguard,1+nguard*k2d:nyb+nguard*k2d, & 
     &            1+nguard*k3d:nzb+nguard*k3d,idest) = & 
     &  unk(ivar,1+nguard0:nxb+nguard0,1+nguard0*k2d:nyb+nguard0*k2d, & 
     &            1+nguard0*k3d:nzb+nguard0*k3d,lb)
              endif
            enddo

          elseif(lb.gt.lnblocks) then

          vtype = 1
          call mpi_set_message_limits( & 
     &                 dtype,ia,ib,ja,jb,ka,kb,vtype)
          index = index0

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
              unk1(ivar_next,i,j,k,idest) = temprecv_buf(index+ivar)
            enddo
            index = index+ngcell_on_cc
          enddo
          enddo
          enddo

          endif

          endif ! end if (lcc

          if(lfc) then

            if(lb.le.lnblocks) then

            do ivar=1,nfacevar
              if(int_gcell_on_fc(1,ivar)) then
            facevarx1(ivar,1+nguard:nxb+nguard+1, & 
     &                  1+nguard*k2d:nyb+nguard*k2d, & 
     &                  1+nguard*k3d:nzb+nguard*k3d,idest) = & 
     &       facevarx(ivar,1+nguard0:nxb+nguard0+1, & 
     &                  1+nguard0*k2d:nyb+nguard0*k2d, & 
     &                  1+nguard0*k3d:nzb+nguard0*k3d,lb)
              endif
            enddo

          elseif(lb.gt.lnblocks) then

! starting index if cell-centered data is also included in recv_buf
          index = index0 + ngcell_on_cc*message_size_cc(dtype)

          vtype = 2
          call mpi_set_message_limits( & 
     &                 dtype,ia,ib,ja,jb,ka,kb,vtype)

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
     &                 temprecv_buf(index+ivar)
            enddo
            index = index+ngcell_on_fc(1)
          enddo
          enddo
          enddo

          endif

            if(ndim.ge.2) then

            if(lb.le.lnblocks) then

            do ivar=1,nfacevar
              if(int_gcell_on_fc(2,ivar)) then
              facevary1(ivar,1+nguard:nxb+nguard, & 
     &                   1+nguard*k2d:nyb+(nguard+1)*k2d, & 
     &                   1+nguard*k3d:nzb+nguard*k3d,idest) = & 
     &        facevary(ivar,1+nguard0:nxb+nguard0, & 
     &                   1+nguard0*k2d:nyb+(nguard0+1)*k2d, & 
     &                   1+nguard0*k3d:nzb+nguard0*k3d,lb)
              endif
            enddo

          elseif(lb.gt.lnblocks) then

          vtype = 3
          call mpi_set_message_limits( & 
     &                 dtype,ia,ib,ja,jb,ka,kb,vtype)

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
     &                 temprecv_buf(index+ivar)
            enddo
            index = index+ngcell_on_fc(2)
          enddo
          enddo
          enddo

          endif

            endif  ! end (if ndim

            if(ndim.eq.3) then

            if(lb.le.lnblocks) then

            do ivar=1,nfacevar
              if(int_gcell_on_fc(3,ivar)) then
              facevarz1(ivar,1+nguard:nxb+nguard, & 
     &                    1+nguard*k2d:nyb+nguard*k2d, & 
     &                    1+nguard*k3d:nzb+(nguard+1)*k3d,idest) = & 
     &         facevarz(ivar,1+nguard0:nxb+nguard0, & 
     &                    1+nguard0*k2d:nyb+nguard0*k2d, & 
     &                    1+nguard0*k3d:nzb+(nguard0+1)*k3d,lb)
              endif
            enddo

            elseif(lb.gt.lnblocks) then

          vtype = 4
          call mpi_set_message_limits( & 
     &                 dtype,ia,ib,ja,jb,ka,kb,vtype)

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
     &                 temprecv_buf(index+ivar)
            enddo
            index = index+ngcell_on_fc(3)
          enddo
          enddo
          enddo

          endif

            endif                 ! end if (ndim

          endif                   ! end of lfc if test

          if (ndim > 1) then
          if(lec) then

            if(lb.le.lnblocks) then

            do ivar=1,nvaredge
              if(int_gcell_on_ec(1,ivar)) then
            unk_e_x1(ivar,1+nguard:nxb+nguard, & 
     &                 1+nguard*k2d:nyb+nguard*k2d+k2d, & 
     &                 1+nguard*k3d:nzb+nguard*k3d+k3d,idest) = & 
     &        unk_e_x(ivar,1+nguard0:nxb+nguard0, & 
     &                  1+nguard0*k2d:nyb+nguard0*k2d+k2d, & 
     &                  1+nguard0*k3d:nzb+nguard0*k3d+k3d,lb)
              endif
            enddo

            elseif(lb.gt.lnblocks) then

! starting index if cell-centered data is also included in recv_buf
         index = index0 + ngcell_on_cc*message_size_cc(dtype) & 
                        + ngcell_on_fc(1) * message_size_fcx(dtype) &
                        + ngcell_on_fc(2) * message_size_fcy(dtype) &
                        + ngcell_on_fc(3) * message_size_fcz(dtype)

          vtype = 5
          call mpi_set_message_limits( & 
     &                 dtype,ia,ib,ja,jb,ka,kb,vtype)

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
              unk_e_x1(ivar_next,i,j,k,idest) = & 
     &                 temprecv_buf(index+ivar)
            enddo
            index = index+ngcell_on_ec(1)
          enddo
          enddo
          enddo

          endif

            if(lb.le.lnblocks) then

            do ivar=1,nvaredge
              if(int_gcell_on_ec(2,ivar)) then
            unk_e_y1(ivar,1+nguard:nxb+nguard+1, & 
     &                 1+nguard*k2d:nyb+nguard*k2d, & 
     &                 1+nguard*k3d:nzb+(nguard+1)*k3d,idest) = & 
     &      unk_e_y(ivar,1+nguard0:nxb+nguard0+1, & 
     &                1+nguard0*k2d:nyb+nguard0*k2d, & 
     &                1+nguard0*k3d:nzb+(nguard0+1)*k3d,lb)
              endif
            enddo

            elseif(lb.gt.lnblocks) then

          vtype = 6
          call mpi_set_message_limits( & 
     &                 dtype,ia,ib,ja,jb,ka,kb,vtype)

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
            do ivar=1,ngcell_on_ec(2)
              ivar_next = gcell_on_ec_pointer(2,ivar)
              unk_e_y1(ivar_next,i,j,k,idest) = & 
     &                 temprecv_buf(index+ivar)
            enddo
            index = index+ngcell_on_ec(2)
          enddo
          enddo
          enddo

          endif

            if (ndim == 3) then

            if(lb.le.lnblocks) then

            do ivar=1,nvaredge
              if(int_gcell_on_ec(3,ivar)) then
            unk_e_z1(ivar,1+nguard:nxb+nguard+1, & 
     &                 1+nguard*k2d:nyb+(nguard+1)*k2d, & 
     &                 1+nguard*k3d:nzb+nguard*k3d,idest) = & 
     &        unk_e_z(ivar,1+nguard0:nxb+nguard0+1, & 
     &                  1+nguard0*k2d:nyb+(nguard0+1)*k2d, & 
     &                  1+nguard0*k3d:nzb+nguard0*k3d,lb)
              endif
            enddo

            elseif(lb.gt.lnblocks) then

          vtype = 7
          call mpi_set_message_limits( & 
     &                 dtype,ia,ib,ja,jb,ka,kb,vtype)

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
            do ivar=1,ngcell_on_ec(3)
              ivar_next = gcell_on_ec_pointer(3,ivar)
              unk_e_z1(ivar_next,i,j,k,idest) = & 
     &                 temprecv_buf(index+ivar)
            enddo
            index = index+ngcell_on_ec(3)
          enddo
          enddo
          enddo

          endif

            end if                ! end if (ndim == 3

          endif                   ! end of lec if test
          end if

          if(lnc) then

            if(lb.le.lnblocks) then

            do ivar=1,nvarcorn
              if(int_gcell_on_nc(ivar)) then
            unk_n1(ivar,1+nguard:nxb+nguard+1, & 
     &               1+nguard*k2d:nyb+(nguard+1)*k2d, & 
     &               1+nguard*k3d:nzb+(nguard+1)*k3d,idest) = & 
     &        unk_n(ivar,1+nguard0:nxb+nguard0+1, & 
     &                1+nguard0*k2d:nyb+(nguard0+1)*k2d, & 
     &                1+nguard0*k3d:nzb+(nguard0+1)*k3d,lb)
              endif
            enddo

            elseif(lb.gt.lnblocks) then

! starting index if cell-centered data is also included in recv_buf
         index = index0 + ngcell_on_cc*message_size_cc(dtype) & 
                        + ngcell_on_fc(1) * message_size_fcx(dtype) &
                        + ngcell_on_fc(2) * message_size_fcy(dtype) &
                        + ngcell_on_fc(3) * message_size_fcz(dtype) & 
     &                  + maxval(ngcell_on_ec(1:ndim))* & 
     &                             message_size_ec(dtype)

          vtype = 8
          call mpi_set_message_limits( & 
     &                 dtype,ia,ib,ja,jb,ka,kb,vtype)

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
            do ivar=1,ngcell_on_nc
              ivar_next = gcell_on_nc_pointer(ivar)
              unk_n1(ivar_next,i,j,k,idest) = & 
     &                 temprecv_buf(index+ivar)
            enddo
            index = index+ngcell_on_nc
          enddo
          enddo
          enddo

          endif

          endif                   ! end of lnc if test

          elseif(iopt.ge.2) then
            iopt0 = iopt-1

            if(lb.le.lnblocks) then

            work1(1+nguard_work:nxb+nguard_work, & 
     &            1+nguard_work*k2d:nyb+nguard_work*k2d, & 
     &            1+nguard_work*k3d:nzb+nguard_work*k3d,idest) = & 
     &  work(1+nguard_work0:nxb+nguard_work0, & 
     &       1+nguard_work0*k2d:nyb+nguard_work0*k2d, & 
     &       1+nguard_work0*k3d:nzb+nguard_work0*k3d,lb,iopt0)

            elseif(lb.gt.lnblocks) then

          vtype = 0
          index = index0
          call mpi_set_message_limits( & 
     &                 dtype,ia,ib,ja,jb,ka,kb,vtype)

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
            work1(i,j,k,idest) =  & 
     &              temprecv_buf(index+1)
            index = index+1
          enddo
          enddo
          enddo

          endif

          endif                 ! end of iopt if test



      return
      end subroutine amr_perm_to_1blk
