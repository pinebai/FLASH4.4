      subroutine amr_prolong_cc(mype,iopt,nlayers)


! $RCSfile: amr_prolong_cc.F90,v $
! $Revision: 1.2 $
! $Date: 2004/09/07 00:58:09 $


!------------------------------------------------------------------------
!
! This routine interpolates data from a parent block to any
! newly created child blocks for data which is defined at cell centers.
!
! This routine calls user supplied routines called prolong_unk_fun
! and prolong_work_fun (depending on which data structure the 
! prolongation operation is to be applied) which contains the detailed 
! definition of the prolongation operator which is to be applied to 
! the data in unk or work.
!
! Written :     Peter MacNeice          January 1997
!------------------------------------------------------------------------

use physicaldata
      use tree
      use workspace
      use paramesh_interfaces
      implicit none
      include  'mpif.h'



      integer mype,iopt,nlayers
      integer isg,ichi,ich
      integer jblock,jchild,jface
      integer ioff,joff,koff

!------------------------------------
! local arrays

      real recv(nvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
      integer ichild(2,mchild)

      integer remote_pe,remote_block

      integer nsend,nrecv
      integer reqr(maxblocks_tr)
      integer statr(MPI_STATUS_SIZE,maxblocks_tr)
      integer reqr2(maxblocks)
      integer statr2(MPI_STATUS_SIZE,maxblocks)
      integer block_int2d,block_int3d,myblockint
      integer istart,jstart,kstart
      integer iend,jend,kend
      integer block_point(maxblocks)
      integer errorcode,ierr

      integer irecv_buf(2,mchild,maxblocks)

      logical newparent(maxblocks)

      save ichild

!------------------------------------


! DEFINE BLOCK INTERIOR 
      if (iopt.eq.1) then
         call MPI_TYPE_VECTOR (nyb+2, & 
     &        nvar*(nxb+2), & 
     &        nvar*iu_bnd, & 
     &        MPI_DOUBLE_PRECISION, & 
     &        block_int2d, & 
     &        ierr)
         myblockint = block_int2d
         if (ndim.eq.3) then
            call MPI_TYPE_HVECTOR (nzb+2, & 
     &           1, & 
     &           nvar*iu_bnd*ju_bnd*8, & 
     &           block_int2d, & 
     &           block_int3d, & 
     &           ierr)
            myblockint = block_int3d
         end if
         istart = nguard
         jstart = nguard*k2d+1-k2d
         kstart = nguard*k3d+1-k3d
         iend   = nguard+nxb+1
         jend   = nguard*k2d+nyb+k2d
         kend   = nguard*k3d+nzb+k3d
      else
         call MPI_TYPE_VECTOR (nyb+2, & 
     &        (nxb+2), & 
     &        iuw, & 
     &        MPI_DOUBLE_PRECISION, & 
     &        block_int2d, & 
     &        ierr)
         myblockint = block_int2d
         if (ndim.eq.3) then
            call MPI_TYPE_HVECTOR (nzb+2, & 
     &           1, & 
     &           iuw*juw*8, & 
     &           block_int2d, & 
     &           block_int3d, & 
     &           ierr)
            myblockint = block_int3d
         end if
         istart = nguard_work
         jstart = nguard_work*k2d+1-k2d
         kstart = nguard_work*k3d+1-k3d
         iend   = nguard_work+nxb+1
         jend   = nguard_work*k2d+nyb+k2d
         kend   = nguard_work*k3d+nzb+k3d
      end if
      if (ndim.eq.3) call MPI_TYPE_COMMIT(block_int2d,ierr)
      call MPI_TYPE_COMMIT(myblockint,ierr)

! make sure that nlayers and iopt are set consistently.
      if(iopt.eq.1.and.nlayers.ne.nguard) then
       if(mype.eq.0) then
          open (unit=30,file='amr_log',status='unknown', & 
     &         position='append')
           write(30,*) 'PARAMESH ERROR !'
           write(30,*) 'Error in prolong - iopt and nlayers'
           write(30,*) 'are not consistent. For iopt=1 you must'
           write(30,*) 'set nlayers=nguard.'
           close(30)
           write(30,*) 'PARAMESH ERROR !'
           write(30,*) 'Error in prolong - iopt and nlayers'
           write(30,*) 'are not consistent. For iopt=1 you must'
           write(30,*) 'set nlayers=nguard.'
       endif
       call Driver_abortFlash("Error in prolong: for iopt=1 set nlayers=nguard")
      elseif(iopt.eq.2.and.nlayers.gt.nguard_work) then
       if(mype.eq.0) then
         open (unit=30,file='amr_log',status='unknown', & 
     &        position='append')
           write(30,*) 'PARAMESH ERROR !'
           write(30,*) 'Error in prolong - iopt and nlayers'
           write(30,*) 'are not consistent. For iopt=2 you must'
           write(30,*) 'set nlayers le nguard_work.'
           close(30)
           write(*,*) 'PARAMESH ERROR !'
           write(*,*) 'Error in prolong - iopt and nlayers'
           write(*,*) 'are not consistent. For iopt=2 you must'
           write(*,*) 'set nlayers le nguard_work.'
       endif
       call Driver_abortFlash("Error in prolong: for iopt=2 set nlayers le nguard_work")
      endif

      newparent(:) = .FALSE.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! COMMUNICATE new_child flags to parent and mark them as newparents

      nrecv = 0
      do isg = 1,lnblocks
         do ichi = 1,nchild
            if (child(1,ichi,isg).gt.0) then
               if (child(2,ichi,isg).ne.mype) then
                  nrecv = nrecv + 1
                  call MPI_IRECV(newparent(isg), & 
     &                 1, & 
     &                 MPI_LOGICAL, & 
     &                 child(2,ichi,isg), & 
     &                 child(1,ichi,isg), & 
     &                 gr_meshComm, & 
     &                 reqr(nrecv), & 
     &                 ierr)
               end if
            end if
         end do
      end do
      
! children send messages to parents if off processor

      nsend = 0
      do isg = 1,lnblocks
         if(parent(1,isg).gt.0) then
            if(parent(2,isg).ne.mype) then
               nsend = nsend + 1
               call MPI_SSEND(newchild(isg), & 
     &              1, & 
     &              MPI_LOGICAL, & 
     &              parent(2,isg), &  ! PE TO SEND TO
     &              isg,        &  ! THIS IS THE TAG
     &              gr_meshComm, & 
     &              ierr)
            end if
         end if
      end do
      
      if (nrecv.gt.0) then
         call MPI_WAITALL (nrecv, reqr, statr, ierr)
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      if (lnblocks.gt.0) then

! COMMUNICATE DATA FROM parents to children

      if (iopt.eq.1) then

! count up no. of receives to post

        nrecv = 0
        do isg = 1,lnblocks
          if (newchild(isg)) then
            if (parent(1,isg).gt.0) then
              if (parent(2,isg).ne.mype) then
                nrecv = nrecv + 1
                call MPI_IRECV( & 
     &               unk(1,istart,jstart,kstart,nrecv+lnblocks), & 
     &               1, & 
     &               myblockint, & 
     &               parent(2,isg), & 
     &               isg, & 
     &               gr_meshComm, & 
     &               reqr(nrecv), & 
     &               ierr)
                call MPI_IRECV(irecv_buf(1,1,nrecv), & 
     &               mchild*2, & 
     &               MPI_INTEGER, & 
     &               parent(2,isg), & 
     &               isg+maxblocks, & 
     &               gr_meshComm, & 
     &               reqr2(nrecv), & 
     &               ierr)
                block_point(isg) = nrecv+lnblocks
              end if
            end if
          end if
        end do
        
! parent sends messages to its children if off processor

        nsend = 0
        do isg = 1,lnblocks
          if (newparent(isg)) then
            do ichi = 1,nchild
              if(child(1,ichi,isg).gt.0) then
                if(child(2,ichi,isg).ne.mype) then
                  nsend = nsend + 1
                  call MPI_SSEND(unk(1,istart,jstart,kstart,isg), & 
     &                 1, & 
     &                 myblockint, & 
     &                 child(2,ichi,isg), &  ! PE TO SEND TO
     &                 child(1,ichi,isg), &  ! THIS IS THE TAG
     &                 gr_meshComm, & 
     &                 ierr)
                  call MPI_SSEND(child(1,1,isg), & 
     &                 mchild*2, & 
     &                 MPI_INTEGER, & 
     &                 child(2,ichi,isg), &  ! PE TO SEND TO
     &                 child(1,ichi,isg)+maxblocks, &  ! THIS IS THE TAG
     &                 gr_meshComm, & 
     &                 ierr)
                end if
              end if
            end do
          end if
        end do

        if (nrecv+lnblocks.gt.maxblocks) then
           open (unit=30,file='amr_log',status='unknown', & 
     &          position='append')
           write (30,*)  & 
     &      ' ERROR in prolong: nrecv+lnblocks > maxblocks '
           close(30)
           print *,' ERROR in prolong: nrecv+lnblocks > maxblocks '
           call Driver_abortFlash("ERROR in prolong: nrecv+lnblocks > maxblocks")
        end if
        
        if (nrecv.gt.0) then
          call MPI_WAITALL (nrecv, reqr, statr, ierr)
          call MPI_WAITALL (nrecv, reqr2, statr2, ierr)
        end if

      else                      ! iopt = 2

        nrecv = 0
        do isg = 1,lnblocks
          if (newchild(isg)) then
            if (parent(1,isg).gt.0) then
              if (parent(2,isg).ne.mype) then
                nrecv = nrecv + 1
                call MPI_IRECV( & 
     &               work(istart,jstart,kstart,nrecv+lnblocks,1), & 
     &               1, & 
     &               myblockint, & 
     &               parent(2,isg), & 
     &               isg, & 
     &               gr_meshComm, & 
     &               reqr(nrecv), & 
     &               ierr)
                call MPI_IRECV(irecv_buf(1,1,nrecv), & 
     &               mchild*2, & 
     &               MPI_INTEGER, & 
     &               parent(2,isg), & 
     &               isg+maxblocks, & 
     &               gr_meshComm, & 
     &               reqr2(nrecv), & 
     &               ierr)
                block_point(isg) = nrecv + lnblocks
              end if
            end if
          end if
        end do
        
! parent sends messages to its children if off processor

        nsend = 0
        do isg = 1,lnblocks
          if (newparent(isg)) then
            do ichi = 1,nchild
              if(child(1,ichi,isg).gt.0) then
                if(child(2,ichi,isg).ne.mype) then
                  nsend = nsend + 1
                  call MPI_SSEND(work(istart,jstart,kstart,isg,1), & 
     &                 1, & 
     &                 myblockint, & 
     &                 child(2,ichi,isg), &  ! PE TO SEND TO
     &                 child(1,ichi,isg), &  ! THIS IS THE TAG
     &                 gr_meshComm, & 
     &                 ierr)
                  call MPI_SSEND(child(1,1,isg), & 
     &                 mchild*2, & 
     &                 MPI_INTEGER, & 
     &                 child(2,ichi,isg), &  ! PE TO SEND TO
     &                 child(1,ichi,isg)+maxblocks, &  ! THIS IS THE TAG
     &                 gr_meshComm, & 
     &                 ierr)
                end if
              end if
            end do
          end if
        end do
        
        if (nrecv+lnblocks.gt.maxblocks) then
           open (unit=30,file='amr_log',status='unknown', & 
     &          position='append')
           write (30,*)  & 
     &      ' ERROR in prolong: nrecv+lnblocks > maxblocks '
           close(30)
           print *,' ERROR in prolong: nrecv+lnblocks > maxblocks '
           call Driver_abortFlash("ERROR in prolong: nrecv+lnblocks > maxblocks")
        end if
        
        if (nrecv.gt.0) then
          call MPI_WAITALL (nrecv, reqr, statr, ierr)
          call MPI_WAITALL (nrecv, reqr2, statr2, ierr)
        end if
       
      end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! cycle through the grid blocks on this processor
      do isg = 1,lnblocks

! Is this a new child block ?
        if(newchild(isg)) then


! Does the parent reside on this processor?
          if(parent(2,isg).eq.mype) then
            jblock = parent(1,isg)
            jchild = 0
            do ich=1,nchild
              if (child(1,ich,jblock).eq.isg.and. & 
     &            child(2,ich,jblock).eq.mype) & 
     &             jchild=ich
            enddo
            if (jchild.eq.0) then
               print *,' ERROR1 !!!!, jchild = 0 ',isg,mype
               do ich = 1,nchild
               print *,child(1,ich,jblock),child(2,ich,jblock)
               end do
            end if
          else
            jblock = block_point(isg) - lnblocks
            jchild = 0
            do ich=1,nchild
              if (irecv_buf(1,ich,jblock).eq.isg.and. & 
     &            irecv_buf(2,ich,jblock).eq.mype) & 
     &             jchild=ich
            enddo
            if (jchild.eq.0) then
               print *,' ERROR2 !!!!, jchild = 0 ',isg,mype
               do ich = 1,nchild
               print *,irecv_buf(1,ich,jblock),irecv_buf(2,ich,jblock)
               end do
            end if
          endif

! identify which child this leaf block represents
          
! compute the offset in the parent block appropriate for this child
          ioff = mod(jchild-1,2)*nxb/2
          joff = mod((jchild-1)/2,2)*nyb/2
          koff = mod((jchild-1)/4,2)*nzb/2

! select complete block for prolongation
!          jface=0

! select block interior only for prolongation
          jface = 7

          recv(:,:,:,:) = 0.
          if(iopt.eq.1) then
! copy its parents data into a temporary buffer array.
            if (parent(2,isg).eq.mype) then
              recv(:,istart:iend,jstart:jend,kstart:kend) =  & 
     &              unk(:,istart:iend,jstart:jend,kstart:kend,jblock)
            else
              recv(:,istart:iend,jstart:jend,kstart:kend) =  & 
     &      unk(:,istart:iend,jstart:jend,kstart:kend,block_point(isg))
            end if

! interpolate data from the parent to the child
            call amr_prolong_unk_fun(recv,isg,ioff,joff,koff, & 
     &           jface,mype)

          elseif(iopt.eq.2) then
! copy its parents data into a temporary buffer array.
            if (parent(2,isg).eq.mype) then
              recv1(istart:iend,jstart:jend,kstart:kend) =  & 
     &         work(istart:iend,jstart:jend,kstart:kend,jblock,1)
            else
              recv1(istart:iend,jstart:jend,kstart:kend) =  & 
     &      work(istart:iend,jstart:jend,kstart:kend,block_point(isg),1)
            end if

! interpolate data from the parent to the child
            call amr_prolong_work_fun(ilw,iuw,jlw,juw, & 
     &           klw,kuw,nlayers, & 
     &           isg,ioff,joff,koff,jface,mype)


          endif

        endif

      enddo

!      end if

      call MPI_TYPE_FREE(myblockint,ierr)
      if (ndim.eq.3) then
         call MPI_TYPE_FREE(block_int2d,ierr)
      end if
      
      return
      end
