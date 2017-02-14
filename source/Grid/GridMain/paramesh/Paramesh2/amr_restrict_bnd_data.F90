      subroutine amr_restrict_bnd_data(mype,idir)


! $RCSfile: amr_restrict_bnd_data.F90,v $
! $Revision: 1.2 $
! $Date: 2004/09/07 00:58:09 $


!------------------------------------------------------------------------
!
! This routine does the data averaging required when a child block
! passes data back to its parent. The parent receives data at the
! block boundary only.
!
! This routine provides a mechanism for passing data defined at block
! boundaries from leaf blocks back to their parents.
! The averaging rules used to combine interface values on the finer
! mesh to construct interface values on the coarser parent mesh are
! specified by the user who provides a function called amr_restrict_red
! to do this.
!
! This routine is only relevant for schemes with even number of grid points.
!
!
! Written :     Peter MacNeice          February 1997
!------------------------------------------------------------------------
use physicaldata
      use tree
      implicit none
      include 'mpif.h'




!------------------------------------
! local arrays

      integer mype,idir

      integer block_point(mchild,maxblocks)
      integer reqr(maxblocks_tr)
      integer statr(MPI_STATUS_SIZE,maxblocks_tr)
      integer remote_pe,remote_block
      integer cnodetype,errorcode

      integer i,j,k,ii,jj,kk,ivar,ierr
      integer istart,iend,jstart,jend,kstart,kend
      integer icoord,len,nrecv,nsend,isg,ichi,ich
      integer jchild
      integer iface,jface,kface,ioff,joff,koff
      
      save cnodetype

!------------------------------------

      if (idir.eq.0) then
        istart = 1
        iend = ndim
      else
        istart = idir
        iend = idir
      end if

      do icoord=istart,iend

      if (icoord.eq.1) then

       len = len_block_bndx*nfluxes
       nrecv = 0
       do isg = 1,lnblocks
         if (nodetype(isg).eq.2) then
           do ichi = 1,nchild
             if (child_type(ichi,isg).eq.1) then
               if(child(2,ichi,isg).ne.mype) then
                 nrecv = nrecv + 1
                 call MPI_IRECV(flux_x(1,1,1,1,lnblocks+nrecv),  & 
     &                  len, & 
     &                  MPI_DOUBLE_PRECISION, & 
     &                  child(2,ichi,isg), & 
     &                  child(1,ichi,isg), & 
     &                  MPI_COMM_WORLD, & 
     &                  reqr(nrecv), & 
     &                  ierr)
                 block_point(ichi,isg) = lnblocks+nrecv
               end if
             end if
           end do
         end if
       end do
       
       if (lnblocks + nrecv .gt. maxblocks) then
         print *,'[AMR_RESTRICT_BND_DATA] ERROR: memory overflow, icoord=1: ', lnblocks, nrecv
	 call Driver_abortFlash("[AMR_RESTRICT_BND_DATA] ERROR: memory overflow, icoord=1")
       end if
       
! children send messages to their parents if off processor

       nsend = 0
       do isg = 1,lnblocks
         if (nodetype(isg).eq.1) then
           if(parent(1,isg).gt.-1) then
             if(parent(2,isg).ne.mype) then
               nsend = nsend + 1
               call MPI_SSEND(flux_x(1,1,1,1,isg), & 
     &              len, & 
     &              MPI_DOUBLE_PRECISION, & 
     &              parent(2,isg), &  ! PE TO SEND TO
     &              isg,        &  ! THIS IS THE TAG
     &              MPI_COMM_WORLD, & 
     &              ierr)
             end if
           end if
         end if
       end do

       elseif (icoord.eq.2) then

       len = len_block_bndy*nfluxes
       nrecv = 0
       do isg = 1,lnblocks
         if (nodetype(isg).eq.2) then
           do ichi = 1,nchild
             if (child_type(ichi,isg).eq.1) then
               if(child(2,ichi,isg).ne.mype) then
                 nrecv = nrecv + 1
                 call MPI_IRECV(flux_y(1,1,1,1,lnblocks+nrecv), & 
     &                  len, & 
     &                  MPI_DOUBLE_PRECISION, & 
     &                  child(2,ichi,isg), & 
     &                  child(1,ichi,isg), & 
     &                  MPI_COMM_WORLD, & 
     &                  reqr(nrecv), & 
     &                  ierr)
                 block_point(ichi,isg) = lnblocks+nrecv
               end if
             end if
           end do
         end if
       end do
      
       if (lnblocks + nrecv .gt. maxblocks) then
         print *,'[AMR_RESTRICT_BND_DATA] ERROR: memory overflow, icoord=2: ', lnblocks, nrecv
	 call Driver_abortFlash("[AMR_RESTRICT_BND_DATA] ERROR: memory overflow, icoord=2")
       end if
      
! children send messages to their parents if off processor

       nsend = 0
       do isg = 1,lnblocks
         if (nodetype(isg).eq.1) then
           if(parent(1,isg).gt.-1) then
             if(parent(2,isg).ne.mype) then
               nsend = nsend + 1
               call MPI_SSEND(flux_y(1,1,1,1,isg), & 
     &              len, & 
     &              MPI_DOUBLE_PRECISION, & 
     &              parent(2,isg), &  ! PE TO SEND TO
     &              isg,        &  ! THIS IS THE TAG
     &              MPI_COMM_WORLD, & 
     &              ierr)
             end if
           end if
         end if
       end do

       elseif (icoord.eq.3) then

       len = len_block_bndz*nfluxes
       nrecv = 0
       do isg = 1,lnblocks
         if (nodetype(isg).eq.2) then
           do ichi = 1,nchild
             if (child_type(ichi,isg).eq.1) then
               if(child(2,ichi,isg).ne.mype) then
                 nrecv = nrecv + 1
                 call MPI_IRECV(flux_z(1,1,1,1,lnblocks+nrecv), & 
     &                  len, & 
     &                  MPI_DOUBLE_PRECISION, & 
     &                  child(2,ichi,isg), & 
     &                  child(1,ichi,isg), & 
     &                  MPI_COMM_WORLD, & 
     &                  reqr(nrecv), & 
     &                  ierr)
                 block_point(ichi,isg) = lnblocks+nrecv
               end if
             end if
           end do
         end if
       end do
       
       if (lnblocks + nrecv .gt. maxblocks) then
         print *,'[AMR_RESTRICT_BND_DATA] ERROR: memory overflow, icoord=3: ', lnblocks, nrecv
	 call Driver_abortFlash("[AMR_RESTRICT_BND_DATA] ERROR: memory overflow, icoord=3")
       end if
      
! children send messages to their parents if off processor

       nsend = 0
       do isg = 1,lnblocks
         if (nodetype(isg).eq.1) then
           if(parent(1,isg).gt.-1) then
             if(parent(2,isg).ne.mype) then
               nsend = nsend + 1
               call MPI_SSEND(flux_z(1,1,1,1,isg), & 
     &              len, & 
     &              MPI_DOUBLE_PRECISION, & 
     &              parent(2,isg), &  ! PE TO SEND TO
     &              isg,        &  ! THIS IS THE TAG
     &              MPI_COMM_WORLD, & 
     &              ierr)
             end if
           end if
         end if
       end do

       end if

       if (nrecv.gt.0) then
         call MPI_WAITALL (nrecv, reqr, statr, ierr)
       end if

! cycle through the grid blocks on this processor
      if(lnblocks.gt.0) then
      do isg = 1,lnblocks

! Is this a parent block of at least one leaf node?
      if(nodetype(isg).eq.2) then

! If yes then cycle through its children.
        do ich=1,nchild

          jchild = ich
          remote_pe = child(2,ich,isg)
          if (remote_pe.eq.mype) then
            remote_block  = child(1,ich,isg)
          else
            remote_block  = block_point(ich,isg)
          end if

! Is this child a leaf block(nodetype=1)?
! If it is then fetch its data.

          cnodetype = child_type(ich,isg)

          if(cnodetype.eq.1) then

              if(icoord.eq.1) then
                 recvarx1(:,:,:,:) = flux_x(:,:,:,:,remote_block)
              elseif(icoord.eq.2) then
                 recvary1(:,:,:,:) = flux_y(:,:,:,:,remote_block)
              elseif(icoord.eq.3) then
                 recvarz1(:,:,:,:) = flux_z(:,:,:,:,remote_block)
              endif


! compute the offset in the parent block appropriate for this child
       iface = mod(jchild-1,2)
       jface = mod((jchild-1)/2,2)
       kface = mod((jchild-1)/4,2)
       ioff = iface*nxb/2
       joff = jface*nyb/2
       koff = kface*nzb/2

! Compute restricted data from the data in the buffer and
! update only boundary values on the parent block
       if(icoord.eq.1) then

         call amr_restrict_red(icoord, isg, ioff, joff, koff)
         i = iface+1
         do k=1+nguard*k3d,nzb+nguard*k3d,2
           kk = (k-nguard*k3d)/2+nguard*k3d+1
           do j=1+nguard*k2d,nyb+nguard*k2d,2
             jj = (j-nguard*k2d)/2+nguard*k2d+1
             do ivar=1,nfluxes
               flux_x(ivar,i,jj+joff,kk+koff,isg) =  & 
     &                                     bndtempx1(ivar,i,j,k)
             enddo
           enddo
         enddo

       elseif(icoord.eq.2) then

         call amr_restrict_red(icoord, isg, ioff, joff, koff)
         j = jface+1
         do k=1+nguard*k3d,nzb+nguard*k3d,2
           kk = (k-nguard*k3d)/2+nguard*k3d+1
           do i=1+nguard,nxb+nguard,2
             ii = (i-nguard)/2+nguard+1
             do ivar=1,nfluxes
               flux_y(ivar,ii+ioff,j,kk+koff,isg)=bndtempy1(ivar,i,j,k)
             enddo
           enddo
         enddo

       elseif(icoord.eq.3) then

         call amr_restrict_red(icoord, isg, ioff, joff, koff)
         k = kface+1
         do j=1+nguard*k2d,nyb+nguard*k2d,2
           jj = (j-nguard*k2d)/2+nguard*k2d+1
           do i=1+nguard,nxb+nguard,2
             ii = (i-nguard)/2+nguard+1
             do ivar=1,nfluxes
               flux_z(ivar,ii+ioff,jj+joff,k,isg)=bndtempz1(ivar,i,j,k)
             enddo
           enddo
         enddo

       endif

      endif

      enddo

      endif

      enddo

      endif

      enddo

      return
      end


