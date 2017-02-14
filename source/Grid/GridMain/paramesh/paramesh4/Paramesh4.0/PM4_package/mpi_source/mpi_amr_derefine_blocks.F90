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

!#define DEBUG

      subroutine amr_derefine_blocks(lnblocks_old,mype)




! By K. Olson (NASA/GSFC and GMU), 11/96

      use paramesh_dimensions
      use physicaldata
      use tree
      use constants
      use paramesh_comm_data

      implicit none

      integer, intent(inout) :: lnblocks_old
      integer, intent(in)    :: mype

      include 'mpif.h'

! local variables --------------------------------------------------------

      integer :: new_loc(maxblocks_tr)
      integer :: i,j,k,jsend
      integer :: lnblocks2
      integer :: neight(2,mfaces,maxblocks_tr)
      integer :: childt(2,mchild,maxblocks_tr)
      integer :: parentt(2,maxblocks_tr)
      integer :: statr(MPI_STATUS_SIZE,maxblocks_tr)
      integer :: reqr(maxblocks_tr)
      integer :: ierr,nsend,nrecv
      integer :: nodetype_chi(nchild,maxblocks_tr)

      logical :: lsend
      integer :: jr0
      real    :: eps,accuracy

! ------------------------------------------------------------------------


#ifdef DEBUG
       Call MPI_BARRIER(amr_mpi_meshComm, ierr)
       write(*,*) 'Entering-derefine-blocks-pe=',mype
       write(*,*) 'Entering derefine blocks pe=',mype
       write(*,*) 'Entering derefine blocks pe=',mype
       Call MPI_BARRIER(amr_mpi_meshComm, ierr)
#endif /* DEBUG */

      accuracy = 10./10.**precision(accuracy)
      eps = accuracy

! remove blocks marked for derefinement by packing the data

      do i = 1,maxblocks_tr
         new_loc(i) = -1
      end do

! Compute new_loc, new_loc marks where each block will end up after the
! derefinement is done

      k = 1
      do i = 1,lnblocks
         if (.not.derefine(i)) then
            new_loc(i) = k
            k = k + 1
          endif
      end do

! 4) reconnect all pointers

      parentt(:,1:lnblocks) = parent(:,1:lnblocks)
      childt(:,:,1:lnblocks) = child(:,:,1:lnblocks)
      neight(:,:,1:lnblocks) = neigh(:,:,1:lnblocks)

      nrecv = 0
      do i = 1,lnblocks
         if (parent(1,i).gt.0) then
           if (parent(2,i).ne.mype) then
             nrecv = nrecv + 1
             call MPI_int_IRECV(parentt(1,i),1,MPI_INTEGER, & 
     &            parent(2,i),i,amr_mpi_meshComm, & 
     &            reqr(nrecv),ierr)
           else
             parentt(1,i) = new_loc(parent(1,i))
           end if
         end if
       end do
       
       nsend = 0
       do i = 1,lnblocks
         do j = 1,nchild
           if (child(1,j,i).gt.0) then
             if (child(2,j,i).ne.mype) then
               ! parent is sending to all its children
               nsend = nsend + 1
               call MPI_int_SSEND (new_loc(i),1,MPI_INTEGER, & 
     &              child(2,j,i),child(1,j,i),amr_mpi_meshComm, & 
     &              ierr)
             end if
           end if
         end do
       end do

      if (nrecv.gt.0) then
        call MPI_WAITALL(nrecv,reqr,statr,ierr)
      end if

      nrecv = 0
      do i = 1,lnblocks
        do j = 1,nchild
          if (child(1,j,i).gt.0) then
            if (child(2,j,i).ne.mype) then
              nrecv = nrecv + 1
              call MPI_int_IRECV(childt(1,j,i),1,MPI_INTEGER, & 
     &             child(2,j,i),child(1,j,i),amr_mpi_meshComm, & 
     &             reqr(nrecv),ierr)
            else
              childt(1,j,i) = new_loc(child(1,j,i))
            end if
          end if
        end do
       end do
       
       nsend = 0
       do i = 1,lnblocks
         if (parent(1,i).gt.0) then
           if (parent(2,i).ne.mype) then
! child is sending to its parent
             nsend = nsend + 1
             call MPI_int_SSEND (new_loc(i),1,MPI_INTEGER, & 
     &            parent(2,i),i,amr_mpi_meshComm, & 
     &            ierr)
           end if
         end if
       end do

      if (nrecv.gt.0) then
        call MPI_WAITALL(nrecv,reqr,statr,ierr)
      end if

      do j = 1,nfaces
#ifdef DEBUG
       Call MPI_BARRIER(amr_mpi_meshComm, ierr)
         write(*,*) 'derefineblocks : pe ',mype,' j ',j
       Call MPI_BARRIER(amr_mpi_meshComm, ierr)
#endif /* DEBUG */

         if (mod(j,2).eq.0) then
            jsend = j - 1
         else
            jsend = j + 1
         end if
            
         nrecv = 0
         do i = 1,lnblocks
            if (neigh(1,j,i).gt.0) then
               if (neigh(2,j,i).ne.mype) then
                  nrecv = nrecv + 1
                  call MPI_int_IRECV(neight(1,j,i),1,MPI_INTEGER, & 
     &                 neigh(2,j,i),neigh(1,j,i),amr_mpi_meshComm, & 
     &                 reqr(nrecv),ierr)
#ifdef DEBUG
         write(*,*) 'derefineblocks : pe ',mype,' blk ',i,' post recv ', & 
     &               nrecv,' tag ',neigh(1,j,i),' src pe ',neigh(2,j,i)
#endif /* DEBUG */
               else
                  neight(1,j,i) = new_loc(neigh(1,j,i))
               end if
            end if
         end do
      
         nsend = 0
         do i = 1,lnblocks

         lsend=.true.

         if(spherical_pm) then
         if(j.eq.3.and.abs(bnd_box(2,2,i)-pi).lt.eps) lsend=.false.
         if(j.eq.4.and.abs(bnd_box(1,2,i)).lt.eps) lsend=.false.
         endif

          if(lsend) then

           if (neigh(1,jsend,i).gt.0) then
             if (neigh(2,jsend,i).ne.mype) then
               nsend = nsend + 1
#ifdef DEBUG
       write(*,*) 'derefineblocks : pe ',mype,' blk ',i,' post n send ', & 
     &               nsend ,' tag ',i,' jsend ',jsend, & 
     &              ' bnd_box ',bnd_box(:,2,i), & 
     &              ' neigh(:,jr0,i) ',neigh(:,jsend,i)
#endif /* DEBUG */
               call MPI_int_SSEND (new_loc(i),1,MPI_INTEGER, & 
     &              neigh(2,jsend,i),i,amr_mpi_meshComm, & 
     &              ierr)
             end if
           end if

         end if

         if(spherical_pm) then
           lsend = .true.
           jr0 = jsend
           if(j.eq.3.and.abs(bnd_box(1,2,i)).lt.eps) jr0 = 3
           if(j.eq.4.and.abs(bnd_box(1,2,i)).lt.eps) lsend=.false.
           if(j.eq.4.and.abs(bnd_box(2,2,i)-pi).lt.eps) jr0 = 4
           if(j.eq.3.and.abs(bnd_box(2,2,i)-pi).lt.eps) lsend=.false.
           if(abs(bnd_box(1,2,i)).lt.eps.and. & 
     &        abs(bnd_box(2,2,i)-pi).lt.eps) then
             write(*,*) 'both poles in blk ',i,abs(bnd_box(1,2,i)), & 
     &         abs(bnd_box(2,2,i)-pi)
             lsend=.true.
           endif

           if(lsend.and.jr0.eq.j) then

             if (neigh(1,jr0,i).gt.0) then
               if (neigh(2,jr0,i).ne.mype) then
                  nsend = nsend + 1
#ifdef DEBUG
         write(*,*) 'derefineblocks : pe ',mype,' blk ',i,' post send ', & 
     &               nsend ,' tag ',i,' jr0 ',jr0, & 
     &              ' bnd_box ',bnd_box(:,2,i), & 
     &              ' neigh(:,jr0,i) ',neigh(:,jr0,i)
#endif /* DEBUG */
                  call MPI_int_SSEND (new_loc(i),1,MPI_INTEGER, & 
     &                 neigh(2,jr0,i),i,amr_mpi_meshComm, & 
     &                 ierr)
               end if
             end if

            end if
         endif

         end do

#ifdef DEBUG
       Call MPI_BARRIER(amr_mpi_meshComm, ierr)
         write(*,*) 'derefineblocks a: pe ',mype,' nrecv nsend ', & 
     &          nrecv, nsend,' face ',j
       Call MPI_BARRIER(amr_mpi_meshComm, ierr)
#endif /* DEBUG */

         if (nrecv.gt.0) then
            call MPI_WAITALL(nrecv,reqr,statr,ierr)
         end if
#ifdef DEBUG
       Call MPI_BARRIER(amr_mpi_meshComm, ierr)
         write(*,*) 'derefineblocks b: pe ',mype
       Call MPI_BARRIER(amr_mpi_meshComm, ierr)
#endif /* DEBUG */

      end do

      do i = 1,lnblocks_old
        if (new_loc(i).ne.i.and.new_loc(i).gt.0) then
          if(nvar.gt.0) unk(:,:,:,:,new_loc(i)) = unk(:,:,:,:,i)
          if (nfacevar.gt.0) then
             facevarx(:,:,:,:,new_loc(i)) = facevarx(:,:,:,:,i)
             facevary(:,:,:,:,new_loc(i)) = facevary(:,:,:,:,i)
             facevarz(:,:,:,:,new_loc(i)) = facevarz(:,:,:,:,i)
          end if
          if (nvaredge.gt.0) then
             unk_e_x(:,:,:,:,new_loc(i)) = unk_e_x(:,:,:,:,i)
             unk_e_y(:,:,:,:,new_loc(i)) = unk_e_y(:,:,:,:,i)
             unk_e_z(:,:,:,:,new_loc(i)) = unk_e_z(:,:,:,:,i)
          end if
          if (nvarcorn.gt.0) unk_n(:,:,:,:,new_loc(i)) =  & 
     &                                       unk_n(:,:,:,:,i)
        end if
      end do

      parent(1,1:lnblocks) = parentt(1,1:lnblocks)
      child(1,:,1:lnblocks) = childt(1,:,1:lnblocks)
      neigh(1,:,1:lnblocks) = neight(1,:,1:lnblocks)

      k = 1
      lnblocks2 = lnblocks
      do i = 1,lnblocks
         
         if (.not.derefine(i)) then
            
            if (k.ne.i) then
               do j = 1,nchild
                  child(1,j,k) = child(1,j,i)
                  child(2,j,k) = child(2,j,i)
               end do
               parent(1,k) = parent(1,i)
               parent(2,k) = parent(2,i)
               do j = 1,nfaces
                  neigh(1,j,k) = neigh(1,j,i)
                  neigh(2,j,k) = neigh(2,j,i)
               end do
               do j = 1,ndim
                  coord(j,k) = coord(j,i)
                  bnd_box(1,j,k) = bnd_box(1,j,i)
                  bnd_box(2,j,k) = bnd_box(2,j,i)
               end do
               bsize(:,k) = bsize(:,i)
               newchild(k) = newchild(i)
               which_child(k) = which_child(i)
               lrefine(k) = lrefine(i)
               bflags(:,k) = bflags(:,i)
               work_block(k) = work_block(i)
               if (empty_cells) then
                  empty(k) = empty(i)
               end if
               
            end if

            k = k + 1
            
         else
            
            lnblocks2 = lnblocks2 - 1
            lnblocks_old = lnblocks_old - 1
            
         end if
         
      end do

! 3) overwrite old locations

      do i = lnblocks2+1,lnblocks
         
         derefine(i) = .FALSE.
         do j = 1,nchild
            child(1,j,i) = -1
            child(2,j,i) = -1
         end do
         parent(1,i) = -1
         parent(2,i) = -1
         do j = 1,nfaces
            neigh(1,j,i) = -1
            neigh(2,j,i) = -1
         end do
         do j = 1,ndim
            coord(j,i) = -1.
            bnd_box(1,j,i) = -1.
            bnd_box(2,j,i) = -1.
         end do
         bsize(:,i) = -1.
         nodetype(i) = -1
         which_child(i) = -1
         newchild(i) = .FALSE.
         lrefine(i) = -1
         bflags(:,i) = -1
         work_block(i) = 0.
         if (empty_cells) then
            empty(i) = 0
         end if
         
      end do
      
      lnblocks = lnblocks2
#ifdef DEBUG
         write(*,*) 'derefineblocks c: pe ',mype
#endif /* DEBUG */


! reset node types

      do i = 1,lnblocks
         nodetype(i) = 3
         if (child(1,1,i).le.-1) then
            nodetype(i) = 1
         end if
      end do
      nrecv = 0
      do i = 1,lnblocks
         do j = 1,nchild
            nodetype_chi(j,i) = -1 
            if (child(1,j,i).gt.-1) then
            if (child(2,j,i).ne.mype) then
               nrecv = nrecv + 1
               call MPI_int_IRECV(nodetype_chi(j,i), & 
     &                        1, & 
     &                        MPI_INTEGER, & 
     &                        child(2,j,i), & 
     &                        child(1,j,i), & 
     &                        amr_mpi_meshComm, & 
     &                        reqr(nrecv), & 
     &                        ierr)
            else
               nodetype_chi(j,i) = nodetype(child(1,j,i))
            end if
            end if
         end do
      end do

      nsend = 0
      do i = 1,lnblocks
! send nodetype to your parent
         if (parent(1,i).ge.1) then
         if (parent(2,i).ne.mype) then
            nsend = nsend + 1
 ! parent(2,i) - PE TO SEND TO
 ! i           - THIS IS THE TAG
            call MPI_int_SSEND(nodetype(i), & 
     &                     1, & 
     &                     MPI_INTEGER, & 
     &                     parent(2,i), & 
     &                     i, & 
     &                     amr_mpi_meshComm, & 
     &                     ierr)
         end if
         end if
      end do

      if (nrecv.gt.0) then
         call MPI_WAITALL (nrecv, reqr, statr, ierr)
      end if

      do i = 1,lnblocks
         do j = 1,nchild
            if (nodetype_chi(j,i).eq.1) nodetype(i) = 2
         end do
      end do

! reset derefine flags

      do i = 1,maxblocks_tr
         derefine(i) = .FALSE.
      end do



      return
      end subroutine amr_derefine_blocks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine amr_check_derefine (mype)

! By K. Olson (NASA/GSFC and GMU), 4/97

      use paramesh_dimensions
      use physicaldata
      use tree
      use timings
      use paramesh_comm_data

      implicit none

      integer, intent(in) :: mype

#ifdef AUTOPACK
      include '/home/macneice/Autopack/1.3.2/include/autopack.fh'
      integer :: count_ap,done,i_ap,tag_ap
      integer :: rbuf, size_ap, sender, rtag, status_ap
      integer :: sdata(2), rdata(2)
#else
      include 'mpif.h'
#endif /* AUTOPACK */

! local variables ----------------------------------------------------------

      integer :: i,j
      integer :: nodetype2(maxblocks_tr)
      integer :: nodetype_recv(maxblocks_tr)
      integer :: tnodetype_recv(maxblocks_tr)
      integer :: nodetype_send(maxblocks_tr)
      integer :: ipar,ipar_proc
      integer :: reqr(maxblocks_tr)
      integer :: statr(MPI_STATUS_SIZE,maxblocks_tr)
      integer :: isg,ierr,neighs,neighr,jsend

      logical :: derefine_chi(nchild,maxblocks_tr)
      logical :: refine_par(maxblocks_tr)
      logical :: trefine_par(maxblocks_tr)

      logical :: unmodified(maxblocks_tr)
      logical :: derefine0(maxblocks_tr)
      integer :: nderefines, loop_count
#ifdef DEBUG
      logical :: derefine_init(maxblocks_tr)
#endif /* DEBUG */

      double precision :: time1

! --------------------------------------------------------------------------
#ifdef DEBUG
      do i = 1,lnblocks
      write(*,*) 'str of check deref: i=',i,' derefine ',derefine(i)
      enddo
      derefine_init = derefine
#endif /* DEBUG */


      no_of_calls_check_derefine = no_of_calls_check_derefine + 1

! If the block is marked for derefinement and it is not
! a leaf block then do not derefine it

      if (timing_mpi) then
         time1 = mpi_wtime()
      endif

      do i = 1,lnblocks
         if (derefine(i).and.nodetype(i).ne.1) derefine(i) = .FALSE.
      end do

#ifdef DEBUG
      do i = 1,lnblocks
        if(derefine_init(i) .and. .not.derefine(i)) & 
     &      write(*,*) 'deref: block i=',i,' switched off at stage 1'
      enddo
      derefine_init = derefine
#endif /* DEBUG */

! Allow the child blocks to only derefine if their parent is NOT
! marked for refinement.
#ifdef AUTOPACK

      call AP_SEND_BEGIN(tag_ap)

      neighs = 0
      do i = 1,lnblocks
         do j = 1,nchild
          if(child(1,j,i).gt.0) then
             if (child(2,j,i).ne.mype) then
                neighs = neighs + 1
! child(2,j,i) - PE TO SEND TO
! child(1,j,i) - THIS IS THE TAG
                sdata(1) = 0
                if(refine(i)) sdata(1) = 1
                sdata(2) = child(1,j,i)
                call AP_BSEND(sdata, & 
     &                        2*4, & 
     &                        child(2,j,i), & 
     &                        child(1,j,i))
!                call MPI_logical_SSEND(refine(i),
!     &                         1,
!     &                         MPI_LOGICAL,
!     &                         child(2,j,i),
!     &                         child(1,j,i),
!     &                         amr_mpi_meshComm,
!     &                         ierr)
            end if
          end if
         end do
      end do

      call AP_FLUSH()
      call AP_SEND_END()


      refine_par(1:lnblocks) = .FALSE.
      i_ap = 0
200   continue
      call AP_RECV_COUNT(count_ap,done)
      if(done.eq.0.or.i_ap.lt.count_ap) then
          call AP_RECV(MPI_ANY_SOURCE,MPI_ANY_TAG, & 
     &                 AP_NOFLAGS, & 
     &                 rbuf, size_ap, sender, rtag, & 
     &                 status_ap)
          if(status_ap.ne.0) then
            call AP_COPY_FREE(rdata(1),rbuf,size_ap)
            refine_par(rdata(2)) = .false.
            if(rdata(1).eq.1) refine_par(rdata(2)) = .true.
            i_ap = i_ap + 1
          endif
          goto 200
      endif

      neighr = i_ap

      call AP_CHECK_SENDS(AP_WAITDEFER,count_ap)
      do i = 1,lnblocks
        if (parent(1,i).gt.0) then
          if (parent(2,i).ne.mype) then
!               neighr = neighr + 1
!               call MPI_logical_IRECV(refine_par(i),
!     .                        1,
!     .                        MPI_LOGICAL,
!     .                        parent(2,i),
!     .                        i,
!     .                        amr_mpi_meshComm,
!     .                        reqr(neighr),
!     .                        ierr)
          else
            refine_par(i) = refine(parent(1,i))
          end if
        end if
      end do

      Call MPI_BARRIER(amr_mpi_meshComm, ierr)

#else /* AUTOPACK */

      neighr = 0
      do i = 1,lnblocks
         refine_par(i) = .FALSE.
         if (parent(1,i).gt.0) then
            if (parent(2,i).ne.mype) then
               neighr = neighr + 1
               call MPI_logical_IRECV(refine_par(i), & 
     &                        1, & 
     &                        MPI_LOGICAL, & 
     &                        parent(2,i), & 
     &                        i, & 
     &                        amr_mpi_meshComm, & 
     &                        reqr(neighr), & 
     &                        ierr)
            else
               refine_par(i) = refine(parent(1,i))
            end if
         end if
      end do

      neighs = 0
      do i = 1,lnblocks
         do j = 1,nchild
          if(child(1,j,i).gt.0) then
             if (child(2,j,i).ne.mype) then
                neighs = neighs + 1
! child(2,j,i) - PE TO SEND TO                
! child(1,j,i) - THIS IS THE TAG
                call MPI_logical_SSEND(refine(i), & 
     &                         1, & 
     &                         MPI_LOGICAL, & 
     &                         child(2,j,i), & 
     &                         child(1,j,i), & 
     &                         amr_mpi_meshComm, & 
     &                         ierr)
            end if
          end if
         end do
      end do

      if (neighr.gt.0) then
         call MPI_WAITALL (neighr, reqr, statr, ierr)
      end if

#endif /* AUTOPACK */


      if (timing_mpi) then
      timer_amr_check_derefine(1) = timer_amr_check_derefine(1) & 
     &                + mpi_wtime() - time1
      time1 = mpi_wtime()
      endif

      do i = 1,lnblocks
         if(nodetype(i).eq.1.and.derefine(i)) then
            if(refine_par(i)) derefine(i)=.false.
         endif
      enddo

#ifdef DEBUG
      do i = 1,lnblocks
        if(derefine_init(i) .and. .not.derefine(i)) & 
     &      write(*,*) 'deref: block i=',i,' switched off at stage 2'
      enddo
      derefine_init = derefine
#endif /* DEBUG */

! Turn off refine flags of non-leaf blocks

      do i = 1,lnblocks
         if (nodetype(i).gt.1.and.refine(i)) refine(i) = .false.
      end do

! Check neighbors to check if OK to derefine

! set nodetype2 = 2 if it either has children or it is marked for
! refinement

! First initialize the record of parents of leaf blocks whose children have all had
! their derefine flags left unmodified.
      unmodified = .true.
      loop_count = 0
      derefine0 = derefine
      nderefines = 1

      do while (loop_count.lt.2.and.nderefines.gt.0)

      if(loop_count.eq.1) derefine = derefine0

      do i = 1,lnblocks
        nodetype2(i) = 1
        if ( (child(1,1,i).ge.1.or.refine(i)) .and. & 
     &        unmodified(i)                ) then ! this node has children 
!        if (child(1,1,i).ge.1.or.refine(i)) then ! this node has children 
                                                 ! or it is marked for 
                                                 ! refinement then its
                                                 ! type is 2
          nodetype2(i) = 2
        end if
      end do

      

! Check for neighboring blocks which are more than one level of refinement
! different


#ifdef AUTOPACK

! cycle through block faces

      do j = 1,nfaces
!      write(*,*) 'pe ',mype,' entered autopack section j= ',j

         if (j.eq.1) jsend = 2
         if (j.eq.2) jsend = 1
         if (j.eq.3) jsend = 4
         if (j.eq.4) jsend = 3
         if (j.eq.5) jsend = 6
         if (j.eq.6) jsend = 5


         call AP_SEND_BEGIN(tag_ap)

! send nodetype2 to neigh if neighbor is off processor and nodetype2 = 2

         neighs = 0
         do isg = 1,lnblocks
!            if (nodetype2(isg).eq.2) then
               if(neigh(1,jsend,isg).gt.-1) then
                  if(neigh(2,jsend,isg).ne.mype) then
                     neighs = neighs + 1
! neigh(2,jsend,isg) - PE TO SEND TO
! isg                - THIS IS THE TAG
                     sdata(1) = nodetype2(isg)
                     sdata(2) = neigh(1,jsend,isg)
                     call AP_BSEND(sdata, & 
     &                              2*4, & 
     &                              neigh(2,jsend,isg), & 
     &                              isg)
!          write(*,*) 'pe ',mype,' sent mess to ',neigh(2,jsend,isg),
!     .               ' with tag ',isg,' data ',nodetype2(isg)
!                  else
!          write(*,*) 'pe ',mype,' no sent mess to ',neigh(2,jsend,isg),
!     .               ' local tag ',isg,' data ',nodetype2(isg)
                  end if
               end if
!            end if
         end do

        call AP_FLUSH()
        call AP_SEND_END()

!       write(*,*) 'pe ',mype,' sends set up autopack section j= ',j

        nodetype_recv(1:lnblocks) = 0
        i_ap = 0
100     continue
        call AP_RECV_COUNT(count_ap,done)
        if(done.eq.0.or.i_ap.lt.count_ap) then
            call AP_RECV(MPI_ANY_SOURCE,MPI_ANY_TAG, & 
     &                   AP_NOFLAGS, & 
     &                   rbuf, size_ap, sender, rtag, & 
     &                   status_ap)
          if(status_ap.ne.0) then
            call AP_COPY_FREE(rdata(1),rbuf,size_ap)
            nodetype_recv(rdata(2)) = rdata(1)
            i_ap = i_ap + 1
          endif
          goto 100
        endif


        call AP_CHECK_SENDS(AP_WAITDEFER,count_ap)

         neighr = i_ap
         mess_counter_chk_deref = mess_counter_chk_deref + i_ap
!        neighr = 0
         do isg = 1,lnblocks
!           nodetype_recv(isg) = 0
            if(neigh(1,j,isg).gt.-1) then
               if(neigh(2,j,isg).ne.mype) then
!                  neighr = neighr + 1
!                  call MPI_int_IRECV(nodetype_recv(isg),
!     .                           1,
!     .                           MPI_INTEGER,
!     .                           neigh(2,j,isg),
!     .                           neigh(1,j,isg),
!     .                           amr_mpi_meshComm,
!     .                           reqr(neighr),
!     .                           ierr)
!                  mess_counter_chk_deref = mess_counter_chk_deref + 1
               else
                  nodetype_recv(isg) = nodetype2(neigh(1,j,isg))
               end if
            end if
         end do

         Call MPI_BARRIER(amr_mpi_meshComm, ierr)

      
        do i = 1,lnblocks
           if (nodetype_recv(i).eq.2) nodetype2(i) = 2
        end do

!      write(*,*) 'pe ',mype,' exited autopack section j=',j
      enddo    ! end loop over faces


#else
! cycle through block faces
      
      do j = 1,nfaces

         if (j.eq.1) jsend = 2
         if (j.eq.2) jsend = 1
         if (j.eq.3) jsend = 4
         if (j.eq.4) jsend = 3
         if (j.eq.5) jsend = 6
         if (j.eq.6) jsend = 5
         
         neighr = 0
         do isg = 1,lnblocks
            nodetype_recv(isg) = 0
            if(neigh(1,j,isg).gt.-1) then
               if(neigh(2,j,isg).ne.mype) then
                  neighr = neighr + 1
                  call MPI_int_IRECV(nodetype_recv(isg), & 
     &                           1, & 
     &                           MPI_INTEGER, & 
     &                           neigh(2,j,isg), & 
     &                           neigh(1,j,isg), & 
     &                           amr_mpi_meshComm, & 
     &                           reqr(neighr), & 
     &                           ierr)
                  mess_counter_chk_deref = mess_counter_chk_deref + 1
               else
                  nodetype_recv(isg) = nodetype2(neigh(1,j,isg))
               end if
            end if
         end do

! send nodetype2 to neigh if neighbor is off processor and nodetype2 = 2

         neighs = 0
         do isg = 1,lnblocks
!            if (nodetype2(isg).eq.2) then
               if(neigh(1,jsend,isg).gt.-1) then
                  if(neigh(2,jsend,isg).ne.mype) then
                     neighs = neighs + 1
! neigh(2,jsend,isg) - PE TO SEND TO
! isg                - THIS IS THE TAG
                     call MPI_int_SSEND(nodetype2(isg), & 
     &                              1, & 
     &                              MPI_INTEGER, & 
     &                              neigh(2,jsend,isg), & 
     &                              isg, & 
     &                              amr_mpi_meshComm, & 
     &                              ierr)
                  end if
               end if
!            end if
         end do
         
         if (neighr.gt.0) then
            call MPI_WAITALL (neighr, reqr, statr, ierr)
         end if
      
         do i = 1,lnblocks
            if (nodetype_recv(i).eq.2) nodetype2(i) = 2
         end do

      end do          ! end of j loop

#endif /* AUTOPACK */


      if (timing_mpi) then
      timer_amr_check_derefine(2) = timer_amr_check_derefine(2) & 
     &                + mpi_wtime() - time1
      time1 = mpi_wtime()
      endif

! Now reset derefine flags based on value of nodetype2

      do i = 1,lnblocks

         if (nodetype2(i).eq.2 .and. derefine(i)) then
            derefine(i) = .FALSE.
         end if
         
      end do

#ifdef DEBUG
      do i = 1,lnblocks
        if(derefine_init(i) .and. .not.derefine(i)) & 
     &      write(*,*) 'deref: block i=',i,' switched off at stage 3', & 
     &                 ' iteration loop_count=',loop_count
      enddo
      derefine_init = derefine
#endif /* DEBUG */

! 1.2) If a block does not have a parent (i.e. = -1) then you can^t derefine
!      it further so if it is marked for derefinement turn derefine off

      do i = 1,lnblocks

         if (derefine(i).and.parent(1,i).lt.0) derefine(i) = .FALSE.

      end do

#ifdef DEBUG
      do i = 1,lnblocks
        if(derefine_init(i) .and. .not.derefine(i)) & 
     &      write(*,*) 'deref: block i=',i,' switched off at stage 4'
      enddo
      derefine_init = derefine
#endif /* DEBUG */

! 1.3) check if all siblings are also marked for derefinement, if not then
!      don^t derefine this block

! parents collect messages from children and count the number of children
! marked for derefinement (stored in nodetype_recv).

      neighr = 0
      do isg = 1,lnblocks
         do j = 1,nchild
            derefine_chi(j,isg) = .FALSE.
            if (child(1,j,isg).gt.-1) then
            if (child(2,j,isg).ne.mype) then
               neighr = neighr + 1
! derefine_chi(j,isg) - this is just junk
               call MPI_logical_IRECV(derefine_chi(j,isg), & 
     &                        1, & 
     &                        MPI_LOGICAL, & 
     &                        child(2,j,isg), & 
     &                        child(1,j,isg), & 
     &                        amr_mpi_meshComm, & 
     &                        reqr(neighr), & 
     &                        ierr)
            else
               derefine_chi(j,isg) = derefine(child(1,j,isg))
            end if
            end if
         end do
      end do

! Children send a message to parent if they are marked for derefinement

      neighs = 0
      nodetype_recv(:) = 0    ! using this variable as a counter here

      do i = 1,lnblocks
!         if (derefine(i)) then
            ipar = parent(1,i) ! parent of i
            ipar_proc = parent(2,i) ! processor parent is stored on
            if (ipar.gt.-1) then
            if (ipar_proc.ne.mype) then
               neighs = neighs + 1
! ipar_proc - PE TO SEND TO
! i         - THIS IS THE TAG
               call MPI_logical_SSEND(derefine(i), & 
     &                        1, & 
     &                        MPI_LOGICAL, & 
     &                        ipar_proc, & 
     &                        i, & 
     &                        amr_mpi_meshComm, & 
     &                        ierr)
            end if
            end if
!         end if
      end do

      if (neighr.gt.0) then
         call MPI_WAITALL (neighr, reqr, statr, ierr)
      end if
      
      if (timing_mpi) then
      timer_amr_check_derefine(3) = timer_amr_check_derefine(3) & 
     &                + mpi_wtime() - time1
      time1 = mpi_wtime()
      endif

      do i = 1,lnblocks
         do j = 1,nchild
            if (derefine_chi(j,i)) then
               nodetype_recv(i) = nodetype_recv(i) + 1
            end if
         end do
      end do
      nodetype_send(1:lnblocks) = nodetype_recv(1:lnblocks)



! Now parent sends nodetype_recv to its children if nodetype_recv = nchild

#ifdef AUTOPACK

      call AP_SEND_BEGIN(tag_ap)

      neighs = 0
      do isg = 1,lnblocks
         if(nodetype(isg).eq.2.and.nodetype_send(isg).eq.nchild & 
     &        .and. unmodified(isg)) then
                      unmodified(isg) = .false.
            endif
            do j = 1,nchild
               if (child(1,j,isg).ge.1) then
                  if (child(2,j,isg).ne.mype) then
                     neighs = neighs + 1
! child(2,j,isg) - PE TO SEND TO
! child(1,j,isg) - THIS IS THE TAG
                sdata(1) = nodetype_send(isg)
                sdata(2) = child(1,j,isg)
                call AP_BSEND(sdata, & 
     &                        2*4, & 
     &                        child(2,j,isg), & 
     &                        child(1,j,isg))
!                     call MPI_int_SSEND(nodetype_send(isg),
!     &                              1,
!     &                              MPI_INTEGER,
!     &                              child(2,j,isg),
!     &                              child(1,j,isg),
!     &                              amr_mpi_meshComm,
!     &                              ierr)
                  else
                     nodetype_recv(child(1,j,isg)) = nodetype_send(isg)
                  end if
               end if
            end do
!         end if
      end do

      call AP_FLUSH()
      call AP_SEND_END()


      i_ap = 0
300   continue
      call AP_RECV_COUNT(count_ap,done)
      if(done.eq.0.or.i_ap.lt.count_ap) then
          call AP_RECV(MPI_ANY_SOURCE,MPI_ANY_TAG, & 
     &                 AP_NOFLAGS, & 
     &                 rbuf, size_ap, sender, rtag, & 
     &                 status_ap)
          if(status_ap.ne.0) then
            call AP_COPY_FREE(rdata(1),rbuf,size_ap)
            nodetype_recv(rdata(2)) = rdata(1)
            i_ap = i_ap + 1
          endif
          goto 300
      endif

      neighr = i_ap

      call AP_CHECK_SENDS(AP_WAITDEFER,count_ap)

      Call MPI_BARRIER(amr_mpi_meshComm, ierr)

#else /* AUTOPACK */

! record modifications to the derefine state of parents, so this
! info can be used in any further loop_count iterations
      do isg = 1,lnblocks
            if(nodetype(isg).eq.2.and.nodetype_send(isg).eq.nchild & 
     &         .and. unmodified(isg)) then
                            unmodified(isg) = .false.
            endif
      enddo

! child blocks post recieves      

      neighr = 0
      do isg = 1,lnblocks
         if(parent(1,isg).gt.-1) then
            if(parent(2,isg).ne.mype) then
               neighr = neighr + 1
               call MPI_int_IRECV(nodetype_recv(isg), & 
     &                        1, & 
     &                        MPI_INTEGER, & 
     &                        parent(2,isg), & 
     &                        isg, & 
     &                        amr_mpi_meshComm, & 
     &                        reqr(neighr), & 
     &                        ierr)
            end if
         end if
      end do

      neighs = 0
      do isg = 1,lnblocks
!         if (nodetype_recv(isg).eq.nchild) then
            do j = 1,nchild
               if (child(1,j,isg).ge.1) then
                  if (child(2,j,isg).ne.mype) then
                     neighs = neighs + 1
! child(2,j,isg) - PE TO SEND TO
! child(1,j,isg) - THIS IS THE TAG
                     call MPI_int_SSEND(nodetype_send(isg), & 
     &                              1, & 
     &                              MPI_INTEGER, & 
     &                              child(2,j,isg), & 
     &                              child(1,j,isg), & 
     &                              amr_mpi_meshComm, & 
     &                              ierr)
                  else
                     nodetype_recv(child(1,j,isg)) = nodetype_send(isg)
                  end if
               end if
            end do
!         end if
      end do
      
      if (neighr.gt.0) then
         call MPI_WAITALL (neighr, reqr, statr, ierr)
      end if

!      do i=1,lnblocks
!        if(nodetype_recv(i).ne.tnodetype_recv(i)) write(*,*) 
!     .       'AUTOERR pe ',mype,' i ',i,' nodetype_recv ',
!     .          nodetype_recv(i),tnodetype_recv(i)
!      enddo

#endif /* AUTOPACK */


! Now loop though the blocks one final time and if nodetype_recv .ne. nchild
!  and
! derefine = .TRUE. then don't derefine

      do isg = 1,lnblocks
         if (derefine(isg).and.nodetype_recv(isg).ne.nchild) then
            derefine(isg) = .FALSE.
         end if
      end do
      
#ifdef DEBUG
      do i = 1,lnblocks
        if(derefine_init(i) .and. .not.derefine(i)) & 
     &      write(*,*) 'deref: block i=',i,' switched off at stage 5'
      enddo
      derefine_init = derefine
#endif /* DEBUG */

      if (timing_mpi) then
      timer_amr_check_derefine(4) = timer_amr_check_derefine(4) & 
     &                + mpi_wtime() - time1
      endif

      nderefines = 0
      do isg = 1,lnblocks
        if( nodetype(isg).eq.1 .and. & 
     &     (derefine0(isg).and.(.not.derefine(isg))) ) & 
     &     nderefines = nderefines + 1
      end do

      call comm_int_sum_to_all(nderefines,nderefines)

      loop_count = loop_count+1

      enddo                      ! end do while

#ifdef DEBUG
      do i = 1,lnblocks
      write(*,*) 'end of check deref: i=',i,' derefine ',derefine(i)
      enddo
#endif /* DEBUG */

      return
      end subroutine amr_check_derefine



