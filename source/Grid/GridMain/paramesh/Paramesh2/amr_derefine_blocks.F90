      subroutine amr_derefine_blocks(lnblocks_old,new_loc2,mype)


! $RCSfile: amr_derefine_blocks.F90,v $
! $Revision: 1.2 $
! $Date: 2004/09/07 00:58:08 $

#define ERRROR_CHECK

! By K. Olson (NASA/GSFC and GMU), 11/96

      use physicaldata
      use tree
      use Grid_data, ONLY : gr_msgbuffer
      implicit none
      include 'mpif.h'



      integer new_loc(maxblocks_tr),new_loct
      integer new_loc2(2,maxblocks_tr)
      integer deref_blocks(maxblocks_tr),nderef
      integer ichi
      integer i,j,k,lb,jsend
      integer mype,lnblocks2
      integer cempty
      integer lnblocks_old
      integer ii,jj,kk,ivar,n_less
      real times,time_exe

      integer neight(2,mfaces,maxblocks_tr)
      integer childt(2,mchild,maxblocks_tr)
      integer parentt(2,maxblocks_tr),lrefinet(maxblocks_tr)
      integer statr(MPI_STATUS_SIZE,maxblocks_tr)
      integer reqr(maxblocks_tr)
      integer ierr,nsend,nrecv
      integer nodetype_chi(nchild,maxblocks_tr)

      integer sprocs(maxblocks_tr), rprocs(maxblocks_tr)
      integer stags (maxblocks_tr), rtags (maxblocks_tr)
      integer svals (maxblocks_tr), rvals (maxblocks_tr)
      logical slvals(maxblocks_tr), rlvals(maxblocks_tr)
      logical flag


! remove blocks marked for derefinement by packing the data

      do i = 1,maxblocks_tr
         new_loc(i) = -1
      end do

!      nderef = 0
!      do lb = 1,lnblocks
!        if (derefine(lb)) then
!          nderef = nderef + 1
!          deref_blocks(nderef) = lb
!       end if
!      end do
!      do i = 1,nderef
!        j = deref_blocks(i)    
!        do lb = 1,lnblocks
!          if (new_loc2(1,lb).gt.new_loc2(1,j)) 
!     &        new_loc2(1,lb) = new_loc2(1,lb) - 1
!        end do
!      end do

! ABOVE IS not yet correct because it does not take into account the fact
! that if a set of blocks derefines, then all the other blocks that will
! be on that blocks processor after redistribution must know that that
! blocks has derefined so that they can readjust their new_locs accordingly.
! The above only looks at blocks on this processor.

! Compute new_loc, new_loc marks where each block will end up after the
! derefinement is done

      k = 1
      do i = 1,lnblocks
         if (.not.derefine(i)) then
            new_loc(i) = k
            k = k + 1
         end if
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
            if (gr_msgbuffer) then
              rprocs(nrecv) = parent(2,i)
              rtags (nrecv) = i
            else 
              call MPI_IRECV(parentt(1,i),1,MPI_INTEGER, & 
     &             parent(2,i),i,MPI_COMM_WORLD, & 
     &             reqr(nrecv),ierr)
            endif
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
               if (gr_msgbuffer) then
                 sprocs(nsend) = child(2,j,i)
                 stags (nsend) = child(1,j,i)
                 svals (nsend) = new_loc(i)
               else
                 call MPI_SSEND (new_loc(i),1,MPI_INTEGER, & 
     &              child(2,j,i),child(1,j,i),MPI_COMM_WORLD, & 
     &              ierr)
               endif
             end if
           end if
         end do
       end do

       if (.not.(gr_msgbuffer)) then
         if (nrecv.gt.0) then
           call MPI_WAITALL(nrecv,reqr,statr,ierr)
         end if
       else
         call b_int_sendrcv( 128,1,nsend,sprocs,stags,svals, & 
     &                             nrecv,rprocs,rtags,rvals)

         nrecv = 0
         do i = 1,lnblocks
           if (parent(1,i).gt.0) then
             if (parent(2,i).ne.mype) then
               nrecv = nrecv + 1
#ifdef ERRROR_CHECK
                    if( (rprocs(nrecv).ne. parent(2,i)).or. & 
     &                  (rtags (nrecv).ne. i) ) then
                       print *,mype,'ASomething horrible has happend.'
                    endif
#endif
               parentt(1,i) = rvals(nrecv)
             end if
           end if
         end do
       endif


      nrecv = 0
      do i = 1,lnblocks
        do j = 1,nchild
          if (child(1,j,i).gt.0) then
            if (child(2,j,i).ne.mype) then
              nrecv = nrecv + 1
              if (gr_msgbuffer) then 
                rprocs(nrecv) = child(2,j,i)
                rtags (nrecv) = child(1,j,i)
              else
                call MPI_IRECV(childt(1,j,i),1,MPI_INTEGER, & 
     &             child(2,j,i),child(1,j,i),MPI_COMM_WORLD, & 
     &             reqr(nrecv),ierr)
              endif
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
             if (gr_msgbuffer) then 
               sprocs(nsend) = parent(2,i) 
               stags (nsend) = i
               svals (nsend) = new_loc(i)
             else
               call MPI_SSEND (new_loc(i),1,MPI_INTEGER, & 
     &            parent(2,i),i,MPI_COMM_WORLD, & 
     &            ierr)
             endif
           end if
         end if
       end do

      if (.not.(gr_msgbuffer)) then
        if (nrecv.gt.0) then
          call MPI_WAITALL(nrecv,reqr,statr,ierr)
        end if
      else

         call b_int_sendrcv( 128,1,nsend,sprocs,stags,svals, & 
     &                             nrecv,rprocs,rtags,rvals)
        nrecv = 0
        do i = 1,lnblocks
          do j = 1,nchild
            if (child(1,j,i).gt.0) then
              if (child(2,j,i).ne.mype) then
                nrecv = nrecv + 1
#ifdef ERRROR_CHECK
                    if( (rprocs(nrecv).ne. child(2,j,i)).or. & 
     &                  (rtags (nrecv).ne. child(1,j,i)) ) then
                       print *,mype,'BSomething horrible has happend.'
                    endif
#endif
                childt(1,j,i) = rvals(nrecv)
              end if
            end if
          end do
        end do
      endif

      nrecv = 0
      nsend = 0
      do j = 1,nfaces

         if (mod(j,2).eq.0) then
            jsend = j - 1
         else
            jsend = j + 1
         end if
            
         if (.not.(gr_msgbuffer)) nrecv = 0

         do i = 1,lnblocks
            if (neigh(1,j,i).gt.0) then
               if (neigh(2,j,i).ne.mype) then
                  nrecv = nrecv + 1
                  if (gr_msgbuffer) then
                    rprocs(nrecv) = neigh(2,j,i)
                    rtags (nrecv) = neigh(1,j,i)
                  else
                    call MPI_IRECV(neight(1,j,i),1,MPI_INTEGER, & 
     &                 neigh(2,j,i),neigh(1,j,i),MPI_COMM_WORLD, & 
     &                 reqr(nrecv),ierr)
                  endif
               else
                  neight(1,j,i) = new_loc(neigh(1,j,i))
               end if
            end if
         end do
      

         if (.not.(gr_msgbuffer)) nsend = 0
         do i = 1,lnblocks
            if (neigh(1,jsend,i).gt.0) then
               if (neigh(2,jsend,i).ne.mype) then
                  nsend = nsend + 1
                  if (gr_msgbuffer) then
                    sprocs(nsend) = neigh(2,jsend,i)
                    stags (nsend) = i
                    svals (nsend) = new_loc(i)
                  else 
                    call MPI_SSEND (new_loc(i),1,MPI_INTEGER, & 
     &                 neigh(2,jsend,i),i,MPI_COMM_WORLD, & 
     &                 ierr)
                  endif
               end if
            end if
         end do

         if (.not.(gr_msgbuffer)) then
           if (nrecv.gt.0) then
              call MPI_WAITALL(nrecv,reqr,statr,ierr)
           end if
         end if

      end do

      if (gr_msgbuffer) then
        call b_int_sendrcv( 128,1,nsend,sprocs,stags,svals, & 
     &                             nrecv,rprocs,rtags,rvals)
        nrecv = 0
        do j = 1,nfaces

           if (mod(j,2).eq.0) then
              jsend = j - 1
           else
              jsend = j + 1
           end if
            
           do i = 1,lnblocks
              if (neigh(1,j,i).gt.0) then
                 if (neigh(2,j,i).ne.mype) then
                    nrecv = nrecv + 1
#ifdef ERRROR_CHECK
                    if( (rprocs(nrecv).ne. neigh(2,j,i)).or. & 
     &                  (rtags (nrecv).ne. neigh(1,j,i)) ) then
                       print *,mype,'CSomething horrible has happend.'
                    endif
#endif
                    neight(1,j,i) = rvals(nrecv)
                  endif
              end if
           end do
        enddo 
      endif



      do i = 1,lnblocks_old
        if (new_loc(i).ne.i.and.new_loc(i).gt.0) then
          unk(:,:,:,:,new_loc(i)) = unk(:,:,:,:,i)
          if (nfacevar.gt.0) then
             facevarx(:,:,:,:,new_loc(i)) = facevarx(:,:,:,:,i)
             facevary(:,:,:,:,new_loc(i)) = facevary(:,:,:,:,i)
             facevarz(:,:,:,:,new_loc(i)) = facevarz(:,:,:,:,i)
          end if
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
!             new_loc2(1,k) = new_loc2(1,i)
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
               lrefine(k) = lrefine(i)
               
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
         newchild(i) = .FALSE.
         lrefine(i) = -1
         
      end do
      
      lnblocks = lnblocks2

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
               if (gr_msgbuffer) then
                 rprocs(nrecv) = child(2,j,i)
                 rtags (nrecv) = child(1,j,i)
               else 
                 call MPI_IRECV(nodetype_chi(j,i), & 
     &                        1, & 
     &                        MPI_INTEGER, & 
     &                        child(2,j,i), & 
     &                        child(1,j,i), & 
     &                        MPI_COMM_WORLD, & 
     &                        reqr(nrecv), & 
     &                        ierr)
               endif
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
            if (gr_msgbuffer) then 
              sprocs(nsend) = parent(2,i)
              stags (nsend) = i
              svals (nsend) = nodetype(i)
            else
              call MPI_SSEND(nodetype(i), & 
     &                     1, & 
     &                     MPI_INTEGER, & 
     &                     parent(2,i), &  ! PE TO SEND TO
     &                     i,       &  ! THIS IS THE TAG
     &                     MPI_COMM_WORLD, & 
     &                     ierr)
            endif
         end if
         end if
      end do

      if (.not.(gr_msgbuffer)) then 
        if (nrecv.gt.0) then
           call MPI_WAITALL (nrecv, reqr, statr, ierr)
        end if
      else
        call b_int_sendrcv( 128,1,nsend,sprocs,stags,svals, & 
     &                            nrecv,rprocs,rtags,rvals)

        nrecv = 0
        do i = 1,lnblocks
           do j = 1,nchild
              if (child(1,j,i).gt.-1) then
              if (child(2,j,i).ne.mype) then
                 nrecv = nrecv + 1
#ifdef ERRROR_CHECK
                    if( (rprocs(nrecv).ne. child(2,j,i)).or. & 
     &                  (rtags (nrecv).ne. child(1,j,i)) ) then
                       print *,mype,'DSomething horrible has happend.'
                    endif
#endif
                 nodetype_chi(j,i) = rvals (nrecv)
              endif
              end if
           end do
        end do
      endif

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
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine amr_check_derefine (mype)

! By K. Olson (NASA/GSFC and GMU), 4/97

      use physicaldata
      use tree
      use Grid_data, ONLY : gr_msgbuffer
      implicit none
      include 'mpif.h'



      integer ichi,ichi_proc
      integer i,j,k,ineigh,ineigh_proc,childt(2),childt2(2,mchild)
      integer cempty
      integer nodetype2(maxblocks_tr)
      integer nodetype_recv(maxblocks_tr)
      integer nodetype_send(maxblocks_tr)
      integer ix,iy,iz,jj,ineigh_old,ineigh_proc_old
      integer ipar,ipar_proc,ineigh2(2)

      integer reqr(maxblocks_tr)
      integer statr(MPI_STATUS_SIZE,maxblocks_tr)
      integer isg,ierr,neighs,neighr,mype,jsend

      logical lt,flag,derefine_chi(nchild,maxblocks_tr)
      logical refine_par(maxblocks_tr)

      integer sprocs(maxblocks_tr), rprocs(maxblocks_tr)
      integer stags (maxblocks_tr), rtags (maxblocks_tr)
      integer svals (maxblocks_tr), rvals (maxblocks_tr)
      logical slvals(maxblocks_tr), rlvals(maxblocks_tr)


! If the block is marked for derefinement and it is not
! a leaf block then do not derefine it

      do i = 1,lnblocks
         if (derefine(i).and.nodetype(i).ne.1) derefine(i) = .FALSE.
      end do

! If a leaf block is marked for derefinement but its parent is
! marked for refinement then cancel the leaf block derefinement.
! DO WE REALLY NEED THIS ????

      neighr = 0
      do i = 1,lnblocks
         refine_par(i) = .FALSE.
         if (parent(1,i).gt.0) then
            if (parent(2,i).ne.mype) then
               neighr = neighr + 1
               if (gr_msgbuffer) then
                 rprocs(neighr) = parent(2,i)
                 rtags (neighr) = i
               else 
                 call MPI_IRECV(refine_par(i), & 
     &                        1, & 
     &                        MPI_LOGICAL, & 
     &                        parent(2,i), & 
     &                        i, & 
     &                        MPI_COMM_WORLD, & 
     &                        reqr(neighr), & 
     &                        ierr)
              endif
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
                if (gr_msgbuffer) then
                  sprocs(neighs) = child(2,j,i)
                  stags (neighs) = child(1,j,i)
                  slvals(neighs) = refine(i)
                else
                  call MPI_SSEND(refine(i), & 
     &                         1, & 
     &                         MPI_LOGICAL, & 
     &                         child(2,j,i),     &  ! PE TO SEND TO
     &                         child(1,j,i),             &  ! THIS IS THE TAG
     &                         MPI_COMM_WORLD, & 
     &                         ierr)
                end if
             end if
          end if
         end do
      end do

      if (.not.(gr_msgbuffer)) then
        if (neighr.gt.0) then
           call MPI_WAITALL (neighr, reqr, statr, ierr)
        endif
      else
        call b_logical_sendrcv(128,1,neighs,sprocs,stags,slvals, & 
     &                               neighr,rprocs,rtags,rlvals)

        neighr = 0
        do i = 1,lnblocks
          if (parent(1,i).gt.0) then
            if (parent(2,i).ne.mype) then
               neighr = neighr + 1
#ifdef ERRROR_CHECK
                    if( (rprocs(neighr).ne. parent(2,i)).or. & 
     &                  (rtags (neighr).ne. i) ) then
                       print *,mype,'ESomething horrible has happend.'
                    endif
#endif
               refine_par(i) = rlvals(neighr)
            end if
          end if
        end do
      end if

      do i = 1,lnblocks
         if(nodetype(i).eq.1.and.derefine(i)) then
            if(refine_par(i)) derefine(i)=.false.
         endif
      enddo

! Check neighbors to check if OK to derefine

! set nodetype2 = 2 if it either has children or it is marked for
! refinement

      do i = 1,lnblocks
        nodetype2(i) = 1
        if (child(1,1,i).ge.1.or.refine(i)) then ! this node has children 
                                                 ! or it is marked for 
                                                 ! refinement then its
                                                 ! type is 2
          nodetype2(i) = 2
        end if
      end do

! Check for neighboring blocks which are more than one level of refinement
! different

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
                  if (gr_msgbuffer) then 
                    rprocs(neighr) = neigh(2,j,isg)
                    rtags (neighr) = neigh(1,j,isg)
                  else
                    call MPI_IRECV(nodetype_recv(isg), & 
     &                           1, & 
     &                           MPI_INTEGER, & 
     &                           neigh(2,j,isg), & 
     &                           neigh(1,j,isg), & 
     &                           MPI_COMM_WORLD, & 
     &                           reqr(neighr), & 
     &                           ierr)
                 endif
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
                     if (gr_msgbuffer) then
                       sprocs(neighs) = neigh(2,jsend,isg)
                       stags (neighs) = isg
                       svals (neighs) = nodetype2(isg)
                     else
                       call MPI_SSEND(nodetype2(isg), & 
     &                              1, & 
     &                              MPI_INTEGER, & 
     &                              neigh(2,jsend,isg), &  ! PE TO SEND TO
     &                              isg,  &  ! THIS IS THE TAG
     &                              MPI_COMM_WORLD, & 
     &                              ierr)
                    end if
                  end if
               end if
!            end if
         end do
         
         if (.not.(gr_msgbuffer)) then
           if (neighr.gt.0) then
             call MPI_WAITALL (neighr, reqr, statr, ierr)
           end if
         else
           call b_int_sendrcv(128,1,neighs,sprocs,stags,svals, & 
     &                              neighr,rprocs,rtags,rvals)
           neighr = 0
           do isg = 1,lnblocks
             if(neigh(1,j,isg).gt.-1) then
                if(neigh(2,j,isg).ne.mype) then
                   neighr = neighr + 1
#ifdef ERRROR_CHECK
                    if( (rprocs(neighr).ne. neigh(2,j,isg)).or. & 
     &                  (rtags (neighr).ne. neigh(1,j,isg)) ) then
                       print *,mype,'FSomething horrible has happend.'
                    endif
#endif
                   nodetype_recv(isg) = rvals(neighr)
                endif
             end if
           end do
         end if

         do i = 1,lnblocks
            if (nodetype_recv(i).eq.2) nodetype2(i) = 2
         end do

      end do

! Now reset derefine flags based on value of nodetype2

      do i = 1,lnblocks

         if (nodetype2(i).eq.2 .and. derefine(i)) then
            derefine(i) = .FALSE.
         end if
         
      end do

! 1.2) If a block does not have a parent (i.e. = -1) then you can^t derefine
!      it further so if it is marked for derefinement turn derefine off

      do i = 1,lnblocks

         if (derefine(i).and.parent(1,i).lt.0) derefine(i) = .FALSE.

      end do

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
               if (gr_msgbuffer) then
                 rprocs(neighr) = child(2,j,isg)    
                 rtags (neighr) = child(1,j,isg)    
               else
                 call MPI_IRECV(derefine_chi(j,isg), &  ! this is just junk
     &                        1, & 
     &                        MPI_LOGICAL, & 
     &                        child(2,j,isg), & 
     &                        child(1,j,isg), & 
     &                        MPI_COMM_WORLD, & 
     &                        reqr(neighr), & 
     &                        ierr)
               endif
            else
               derefine_chi(j,isg) = derefine(child(1,j,isg))
            end if
            end if
         end do
      end do

! Children send a message to parent if they ar marked for derefinement

      neighs = 0
      nodetype_recv(:) = 0    ! using this variable as a counter here

      do i = 1,lnblocks
!         if (derefine(i)) then
            ipar = parent(1,i) ! parent of i
            ipar_proc = parent(2,i) ! processor parent is stored on
            if (ipar.gt.-1) then
            if (ipar_proc.ne.mype) then
               neighs = neighs + 1
               if (gr_msgbuffer) then
                 sprocs(neighs) = ipar_proc
                 stags (neighs) = i
                 slvals(neighs) = derefine(i)
               else
                 call MPI_SSEND(derefine(i), & 
     &                        1, & 
     &                        MPI_LOGICAL, & 
     &                        ipar_proc, &  ! PE TO SEND TO
     &                        i,      &  ! THIS IS THE TAG
     &                        MPI_COMM_WORLD, & 
     &                        ierr)
               endif
            end if
            end if
!         end if
      end do

      if (.not.(gr_msgbuffer)) then
        if (neighr.gt.0) then
           call MPI_WAITALL (neighr, reqr, statr, ierr)
        end if
      else
        call b_logical_sendrcv(128,1,neighs,sprocs,stags,slvals, & 
     &                               neighr,rprocs,rtags,rlvals)
        neighr = 0
        do isg = 1,lnblocks
           do j = 1,nchild
             if (child(1,j,isg).gt.-1) then
             if (child(2,j,isg).ne.mype) then
                neighr = neighr + 1
#ifdef ERRROR_CHECK
                    if( (rprocs(neighr).ne. child(2,j,isg)).or. & 
     &                  (rtags (neighr).ne. child(1,j,isg)) ) then
                       print *,mype,'GSomething horrible has happend.'
                    endif
#endif
                derefine_chi(j,isg) = rlvals(neighr)
            endif
            end if
         end do
        end do
      end if

      do i = 1,lnblocks
         do j = 1,nchild
            if (derefine_chi(j,i)) then
               nodetype_recv(i) = nodetype_recv(i) + 1
            end if
         end do
      end do
      nodetype_send(1:lnblocks) = nodetype_recv(1:lnblocks)

! Now parent sends nodetype_recv to its children if nodetype_recv = nchild

! child blocks post recieves      

      neighr = 0
      do isg = 1,lnblocks
         if(parent(1,isg).gt.-1) then
            if(parent(2,isg).ne.mype) then
               neighr = neighr + 1
               if (gr_msgbuffer) then
                 rprocs(neighr) = parent(2,isg)
                 rtags (neighr) = isg
               else
                 call MPI_IRECV(nodetype_recv(isg), & 
     &                        1, & 
     &                        MPI_INTEGER, & 
     &                        parent(2,isg), & 
     &                        isg, & 
     &                        MPI_COMM_WORLD, & 
     &                        reqr(neighr), & 
     &                        ierr)
              end if
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
                     if (gr_msgbuffer) then
                       sprocs(neighs) = child(2,j,isg)
                       stags (neighs) = child(1,j,isg)
                       svals (neighs) = nodetype_send(isg)
                     else
                       call MPI_SSEND(nodetype_send(isg), & 
     &                              1, & 
     &                              MPI_INTEGER, & 
     &                              child(2,j,isg), &  ! PE TO SEND TO
     &                              child(1,j,isg), &  ! THIS IS THE TAG
     &                              MPI_COMM_WORLD, & 
     &                              ierr)
                    end if
                  else
                     nodetype_recv(child(1,j,isg)) = nodetype_send(isg)
                  end if
               end if
            end do
!         end if
      end do
      
      if (.not.(gr_msgbuffer)) then
        if (neighr.gt.0) then
           call MPI_WAITALL (neighr, reqr, statr, ierr)
        end if
      else
        call b_int_sendrcv(128,1,neighs,sprocs,stags,svals, & 
     &                           neighr,rprocs,rtags,rvals)
        neighr = 0
        do isg = 1,lnblocks
           if(parent(1,isg).gt.-1) then
              if(parent(2,isg).ne.mype) then
                 neighr = neighr + 1
                 nodetype_recv(isg) = rvals(neighr)
              end if
            end if
        end do
      end if

! Now loop though the blocks one final time and if nodetype_recv .ne. nchild
!  and
! derefine = .TRUE. then don't derefine

      do isg = 1,lnblocks
         if (derefine(isg).and.nodetype_recv(isg).ne.nchild) then
            derefine(isg) = .FALSE.
         end if
      end do

      return
      end



!      subroutine neigh_check (nprocs,mype)

!#include "physicaldata.fh"
!      include 'tree.fh'
!      include 'mpif.h'

!      integer neight(2,mfaces)
!      save neight

!      if (mype.eq.0) then
!         open (unit=30,file='amr_log',status='unknown',
!     &        position='append')
!         write (30,*) ' STARTING NEIGHBOR CHECK ',mype
!         close (30)
!      end if

!      do j = 1,nfaces

!         if (mod(j,2).eq.0) then
!            jsend = j - 1
!         else
!            jsend = j + 1
!         end if
!         
!         do i = 1,lnblocks
!            if (neigh(1,j,i).gt.0) then
!               call shmem_integer_get (neight(1,1),
!     &              neigh(1,1,neigh(1,j,i)),2*mfaces,neigh(2,j,i))
!               if (i.eq.8.and.mype.eq.204) then
!               if (neight(1,jsend).ne.i.or.neight(2,jsend).ne.mype) then
!                  open (unit=30,file='amr_log',status='unknown',
!     &                 position='append')
!                  write (30,*) ' ERROR in neighbor: '
!                  write (30,*) ' j = ',j
!                  write (30,*) ' i,mype = ',i,mype
!                  write (30,*) ' neight = ',
!     &                 neight(1,jsend),neight(2,jsend)
!                  write (30,*) ' newchild = ',newchild(i)
!                  do kk = 1,nfaces
!                     write (30,*) ' neigh of i = ',
!     &                    neigh(1,kk,i),neigh(2,kk,i)
!                  end do
!                  do kk = 1,nfaces
!                     write (30,*) ' neigh k = ',
!     &                    neight(1,kk),neight(2,kk)
!                  end do
!                  close(30)
!               end if
!               end if
!            end if
!         end do

!      end do

!      call shmem_barrier_all()
!      if (mype.eq.0) then
!         open (unit=30,file='amr_log',status='unknown',
!     &        position='append')
!         write (30,*) ' DONE NEIGHBOR CHECK ',mype
!         close(30)
!      end if

!      return
!      end


!      subroutine old_loc_check (nprocs,mype,new_loc,old_loc,
!     &     lnblocks2)

!#include "physicaldata.fh"
!      include 'tree.fh'
!      include 'mpif.h'

!      integer neight(2)
!      integer new_loc(2,maxblocks_tr),old_loc(2,maxblocks_tr)
!      save neight

!      if (mype.eq.0) then
!         open (unit=30,file='amr_log',status='unknown',
!     &        position='append')
!         write (30,*) ' STARTING old_loc CHECK ',mype
!         close (30)
!      end if

!      do i = 1,lnblocks
!         if (new_loc(1,i).gt.0) then
!            call shmem_integer_get (neight(1),
!     &              old_loc(1,new_loc(1,i)),2,new_loc(2,i))
!            if (neight(1).ne.i.or.neight(2).ne.mype) then
!               open (unit=30,file='amr_log',status='unknown',
!     &              position='append')
!               write (30,*) ' ERROR in old_loc: '
!               write (30,*) ' i,mype = ',i,mype
!               write (30,*) ' neight = ',
!     &              neight(1),neight(2)
!               close(30)
!            end if
!         end if
!      end do
!      do i = 1,lnblocks2
!         if (old_loc(1,i).gt.0) then
!            call shmem_integer_get (neight(1),
!     &              new_loc(1,old_loc(1,i)),2,old_loc(2,i))
!            if (neight(1).ne.i.or.neight(2).ne.mype) then
!               open (unit=30,file='amr_log',status='unknown',
!     &              position='append')
!               write (30,*) ' ERROR in new_loc: '
!               write (30,*) ' i,mype = ',i,mype
!               write (30,*) ' neight = ',
!     &              neight(1),neight(2)
!               close(30)
!            end if
!         end if
!      end do

!      call shmem_barrier_all()
!      if (mype.eq.0) then
!         open (unit=30,file='amr_log',status='unknown',
!     &        position='append')
!         write (30,*) ' DONE old_loc CHECK ',mype
!         close(30)
!      end if

!      return
!      end
