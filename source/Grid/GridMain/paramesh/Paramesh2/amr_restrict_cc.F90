      subroutine amr_restrict_cc(mype,iopt,rflag)


! $RCSfile: amr_restrict_cc.F90,v $
! $Revision: 1.2 $
! $Date: 2004/09/07 00:58:09 $


!------------------------------------------------------------------------
!
! This routine does the data averaging required when a child block
! passes data back to its parent. The parent receives interior data
! only, not guard cell data.
! This routine calls a user provided routine called restrict_fun
! which defines the pattern of restriction which the user wishes to
! apply.
!
! Written :     Peter MacNeice          December 1996
!------------------------------------------------------------------------

use physicaldata
      use tree
      use workspace
      implicit none
      include 'mpif.h'



      logical rflag(maxblocks)

!------------------------------------
! local arrays
      real send(nvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
      real recv(nvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
      real temp(nvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
      integer block_point(mchild,maxblocks)
      integer reqr(maxblocks_tr)
      integer statr(MPI_STATUS_SIZE,maxblocks_tr)
      integer reqr2(maxblocks)
      integer statr2(MPI_STATUS_SIZE,maxblocks)
      integer errorcode,ierr

      integer mype,iopt
      integer i,j,k,ii,jj,kk,ivar
      integer istart,jstart,kstart
      integer ioff,joff,koff
      integer isg,ich,ichi,neighr,neighs
      integer jchild,jblock
      
      integer remote_pe,remote_block
      integer cnodetype,cempty
      integer block_int2d, block_int3d, myblockint
      save cnodetype,cempty,recv

!------------------------------------

! DEFINE BLOCK INTERIOR 
      call MPI_TYPE_VECTOR (nyb, & 
     &     nvar*nxb, & 
     &     nvar*iu_bnd, & 
     &     MPI_DOUBLE_PRECISION, & 
     &     block_int2d, & 
     &     ierr)
      myblockint = block_int2d
      if (ndim.eq.3) then
        call MPI_TYPE_HVECTOR (nzb, & 
     &       1, & 
     &       nvar*iu_bnd*ju_bnd*8, & 
     &       block_int2d, & 
     &       block_int3d, & 
     &       ierr)
        myblockint = block_int3d
      end if
      istart = nguard+1
      jstart = nguard*k2d+1
      kstart = nguard*k3d+1
      if (ndim.eq.3) call MPI_TYPE_COMMIT(block_int2d,ierr)
      call MPI_TYPE_COMMIT(myblockint,ierr)

      if (lnblocks.gt.0) then

! cycle through the grid blocks on this processor
      do isg = 1,lnblocks

! Is this a leaf block ?
      if(nodetype(isg).eq.1.and.rflag(isg) & 
     &       .and.empty(isg).eq.0) then

! Does the parent reside on this processor?
      if(parent(2,isg).eq.mype) then

! identify which child this leaf block represents
       jchild = 0
       do ich=1,nchild
         if( (child(1,ich,parent(1,isg)).eq.isg) .and. & 
     &       (child(2,ich,parent(1,isg)).eq.mype) )  & 
     &       jchild=ich
       enddo

! compute the offset in the parent block appropriate for this child
       ioff = mod(jchild-1,2)*nxb/2
       joff = mod((jchild-1)/2,2)*nyb/2
       koff = mod((jchild-1)/4,2)*nzb/2

       jblock = parent(1,isg)


! compute its restricted data needed by its parent.
       if(iopt.eq.1) then


       if(mod(nxb,2).eq.0) then
         call amr_restrict_unk_fun(unk(:,:,:,:,isg),temp, jblock, ioff, joff, koff)
       else
         call amr_restrict_unk_fun_recip(unk(:,:,:,:,isg),temp)
       endif


       do k=1+nguard*k3d,nzb+nguard*k3d,2
         kk = (k-nguard*k3d)/2+1+nguard*k3d
         do j=1+nguard*k2d,nyb+nguard*k2d,2
           jj = (j-nguard*k2d)/2+1+nguard*k2d
           do i=1+nguard,nxb+nguard,2
             ii = (i-nguard)/2+1+nguard
             do ivar=1,nvar
               send(ivar,ii,jj,kk) = temp(ivar,i,j,k)
             enddo
           enddo
         enddo
       enddo

       do k=1+nguard*k3d,nzb+(nguard-nzb/2)*k3d
         do j=1+nguard*k2d,nyb+(nguard-nyb/2)*k2d
           do i=1+nguard,nxb+nguard-nxb/2
             do ivar=1,nvar
               unk(ivar,i+ioff,j+joff,k+koff,jblock) = send(ivar,i,j,k)
             enddo
           enddo
         enddo
       enddo

       elseif(iopt.eq.2) then

       if(mod(nxb,2).eq.0) then
         call amr_restrict_work_fun(work(:,:,:,isg,1),temp1, jblock, ioff, joff, koff)
       else
         call amr_restrict_work_fun_recip(work(:,:,:,isg,1),temp1)
       endif

       do k=1+nguard_work*k3d,nzb+nguard_work*k3d,2
         kk = (k-nguard_work*k3d)/2+1+nguard_work*k3d
         do j=1+nguard_work*k2d,nyb+nguard_work*k2d,2
           jj = (j-nguard_work*k2d)/2+1+nguard_work*k2d
           do i=1+nguard_work,nxb+nguard_work,2
             ii = (i-nguard_work)/2+1+nguard_work
             send1(ii,jj,kk) = temp1(i,j,k)
           enddo
         enddo
       enddo

       do k=1+nguard_work*k3d,nzb+(nguard_work-nzb/2)*k3d
         do j=1+nguard_work*k2d,nyb+(nguard_work-nyb/2)*k2d
           do i=1+nguard_work,nxb+nguard_work-nxb/2
             work(i+ioff,j+joff,k+koff,jblock,1) = send1(i,j,k)
           enddo
         enddo
       enddo
       endif

      endif

      endif
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (iopt.eq.1) then

        neighr = 0
        do isg = 1,lnblocks
          if (nodetype(isg).eq.2) then
          do ichi = 1,nchild
            if (child_type(ichi,isg).eq.1) then
              if(child(2,ichi,isg).ne.mype) then
                neighr = neighr + 1
                call MPI_IRECV(unk(1, & 
     &               istart,jstart,kstart,lnblocks+neighr), & 
     &               1, & 
     &               myblockint, & 
     &               child(2,ichi,isg), & 
     &               child(1,ichi,isg), & 
     &               MPI_COMM_WORLD, & 
     &               reqr(neighr), & 
     &               ierr)
                block_point(ichi,isg) = lnblocks+neighr
              end if
            end if
          end do
          end if
        end do
        
        if (lnblocks + neighr .gt. maxblocks) then
           print *,'[AMR_RESTRICT_CC] ERROR: memory overflow, iopt=1: ', neighr, lnblocks
           call Driver_abortFlash("[AMR_RESTRICT_CC] ERROR: memory overflow, iopt=1")
        end if
        
! children send messages to their parents if off processor

        neighs = 0
        do isg = 1,lnblocks
          if (nodetype(isg).eq.1) then
          if(parent(1,isg).gt.-1) then
            if(parent(2,isg).ne.mype) then
              neighs = neighs + 1
              call MPI_SSEND(unk(1,istart,jstart,kstart,isg), & 
     &             1, & 
     &             myblockint, & 
     &             parent(2,isg), &  ! PE TO SEND TO
     &             isg,         &  ! THIS IS THE TAG
     &             MPI_COMM_WORLD, & 
     &             ierr)
            end if
          end if
          end if
        end do
        
        if (neighr.gt.0) then
          call MPI_WAITALL (neighr, reqr, statr, ierr)
        end if
       
      else                      ! iopt = 2

        neighr = 0
        do isg = 1,lnblocks
          if (nodetype(isg).eq.2) then
          do ichi = 1,nchild
            if (child_type(ichi,isg).eq.1) then
              if(child(2,ichi,isg).ne.mype) then
                neighr = neighr + 1
                call MPI_IRECV(work(1,1,1,lnblocks+neighr,1), & 
     &               len_wblock, & 
     &               MPI_DOUBLE_PRECISION, & 
     &               child(2,ichi,isg), & 
     &               child(1,ichi,isg), & 
     &               MPI_COMM_WORLD, & 
     &               reqr(neighr), & 
     &               ierr)
                block_point(ichi,isg) = lnblocks+neighr
               end if
             end if
         end do
         end if
       end do
      
       if (lnblocks + neighr .gt. maxblocks) then
          print *,'[AMR_RESTRICT_CC] ERROR: memory overflow, iopt=2: ', neighr, lnblocks
          call Driver_abortFlash("[AMR_RESTRICT_CC] ERROR: memory overflow, iopt=2")
       end if
      
! children send messages to their parents if off processor

       neighs = 0
       do isg = 1,lnblocks
         if (nodetype(isg).eq.1) then
          if(parent(1,isg).gt.-1) then
           if(parent(2,isg).ne.mype) then
             neighs = neighs + 1
             call MPI_SSEND(work(1,1,1,isg,1), & 
     &            len_wblock, & 
     &            MPI_DOUBLE_PRECISION, & 
     &            parent(2,isg), &  ! PE TO SEND TO
     &            isg,          &  ! THIS IS THE TAG
     &            MPI_COMM_WORLD, & 
     &            ierr)
           end if
         end if
         end if
       end do
       
       if (neighr.gt.0) then
         call MPI_WAITALL (neighr, reqr, statr, ierr)
       end if

      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now parents of leaf nodes which are on remote processors get data
! from their children and then perform restriction on it.

      do isg = 1,lnblocks

! Is this a parent block of at least one leaf node?
      if(nodetype(isg).eq.2.and.rflag(isg)) then

! If yes then cycle through its children.
      do ich=1,nchild

        jchild = ich ! this assumes morton ordering

! Does the child reside on another processor?
        if(child(2,ich,isg).ne.mype) then


! Is this child a leaf block? If it is then fetch its data.
        remote_pe    = child(2,ich,isg)
        remote_block = child(1,ich,isg)
        cnodetype    = child_type(ich,isg)

! NOT USING empty flags for flash code !!!!!

        if(cnodetype.eq.1) then

! compute the offset in the parent block appropriate for this child
! This assumes morton ordering !!!!!
         ioff = mod(jchild-1,2)*nxb/2
         joff = mod((jchild-1)/2,2)*nyb/2
         koff = mod((jchild-1)/4,2)*nzb/2

         if(iopt.eq.1) then

           recv(:,:,:,:) = unk(:,:,:,:,block_point(ich,isg))

! Compute restricted data from the data in the buffer

           jblock = block_point(ich,isg)
           if(mod(nxb,2).eq.0) then
             call amr_restrict_unk_fun(recv,temp, isg, ioff, joff, koff)
           else
             call amr_restrict_unk_fun_recip(recv,temp)
           endif

           do k=1+nguard*k3d,nzb+nguard*k3d,2
             kk = (k-nguard*k3d)/2+1+nguard*k3d
             do j=1+nguard*k2d,nyb+nguard*k2d,2
               jj = (j-nguard*k2d)/2+1+nguard*k2d
               do i=1+nguard,nxb+nguard,2
                 ii = (i-nguard)/2+1+nguard
                 do ivar=1,nvar
                   send(ivar,ii,jj,kk) = temp(ivar,i,j,k)
                 enddo
               enddo
             enddo
           enddo

! update the parent block
           do k=1+nguard*k3d,nzb+(nguard-nzb/2)*k3d
             kk = k+koff
             do j=1+nguard*k2d,nyb+(nguard-nyb/2)*k2d
               jj = j+joff
               do i=1+nguard,nxb+nguard-nxb/2
                 ii = i+ioff
                 do ivar=1,nvar
                   unk(ivar,ii,jj,kk,isg) = & 
     &                          send(ivar,i,j,k)
                 enddo
               enddo
             enddo
           enddo


         elseif(iopt.eq.2) then

           recv1(:,:,:) = work(:,:,:,block_point(ich,isg),1)

! Compute restricted data from the data in the buffer
           if(mod(nxb,2).eq.0) then
             call amr_restrict_work_fun(recv1,temp1, isg, ioff, joff, koff)
           else
             call amr_restrict_work_fun_recip(recv1,temp1)
           endif


           do k=1+nguard_work*k3d,nzb+nguard_work*k3d,2
             kk = (k-nguard_work*k3d)/2+1+nguard_work*k3d
             do j=1+nguard_work*k2d,nyb+nguard_work*k2d,2
               jj = (j-nguard_work*k2d)/2+1+nguard_work*k2d
               do i=1+nguard_work,nxb+nguard_work,2
                 ii = (i-nguard_work)/2+1+nguard_work
                 send1(ii,jj,kk) = temp1(i,j,k)
               enddo
             enddo
           enddo


! update the parent block
           do k=1+nguard_work*k3d,nzb+(nguard_work-nzb/2)*k3d
             kk = k+koff
             do j=1+nguard_work*k2d,nyb+(nguard_work-nyb/2)*k2d
               jj = j+joff
               do i=1+nguard_work,nxb+nguard_work-nxb/2
                 ii = i+ioff
                 work(i+ioff,j+joff,k+koff,isg,1) = send1(i,j,k)
               enddo
             enddo
           enddo

       endif

       endif

       endif

      enddo
      endif
      enddo



      endif

      call MPI_TYPE_FREE(myblockint,ierr)
      if (ndim.eq.3) then
              call MPI_TYPE_FREE(block_int2d,ierr)
      endif

      return
      end


