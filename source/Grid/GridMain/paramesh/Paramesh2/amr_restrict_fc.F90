      subroutine amr_restrict_fc(mype,rflag)


! $RCSfile: amr_restrict_fc.F90,v $
! $Revision: 1.1 $
! $Date: 2004/05/20 01:03:44 $


!------------------------------------------------------------------------
!
! This routine does the data averaging required on the cell face-centered
! data when a child block passes data back to its parent. 
! The parent receives interior data only, not guard cell data.
! This routine calls a user provided routine called restrict_fc_fun
! which defines the pattern of restriction which the user wishes to
! apply.
!
! Written :     Peter MacNeice          July 1997
!------------------------------------------------------------------------

use physicaldata
      use tree
      use workspace
      implicit none


      logical rflag(maxblocks)

!------------------------------------
! local arrays
      real send(nbndvar,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1, & 
     &       kl_bnd:ku_bnd+1)
      real recv(nbndvar,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1, & 
     &       kl_bnd:ku_bnd+1)
      real temp(nbndvar,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1, & 
     &       kl_bnd:ku_bnd+1)

      integer mype
      integer isg,ich,jchild,jblock
      integer ioff,joff,koff
      integer i,j,k,ii,jj,kk,ivar

      integer remote_pe,remote_block
      integer cnodetype,cempty
      save cnodetype,cempty,recv

!------------------------------------



      if(lnblocks.gt.0) then




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
         if( (child(1,ich,parent(1,isg)).eq.isg) .and.  & 
     &       (child(2,ich,parent(1,isg)).eq.mype) ) jchild=ich
       enddo

! compute the offset in the parent block appropriate for this child
       ioff = mod(jchild-1,2)*nxb/2
       joff = mod((jchild-1)/2,2)*nyb/2
       koff = mod((jchild-1)/4,2)*nzb/2

       jblock = parent(1,isg)


! compute its restricted data needed by its parent.

! x face first - note that the indexing of the x face here has a different
! range than for y and z and that the x positions are cell faces while
! the y and z positions are cell centers. These differences are permuted 
! appropriately when we consider the y and z face variables below.

       recv(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd) =  & 
     &  facevarx(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd,isg)
       call amr_restrict_fc_fun(recv,temp,1)

       do k=1+nguard*k3d,nzb+nguard*k3d,2
         kk = (k-nguard*k3d)/2+1+nguard*k3d
         do j=1+nguard*k2d,nyb+nguard*k2d,2
           jj = (j-nguard*k2d)/2+1+nguard*k2d
           do i=1+nguard,nxb+nguard+1,2
             ii = (i-nguard)/2+1+nguard
             do ivar=1,nbndvar
               send(ivar,ii,jj,kk) = temp(ivar,i,j,k)
             enddo
           enddo
         enddo
       enddo

       do k=1+nguard*k3d,nzb+(nguard-nzb/2)*k3d
         do j=1+nguard*k2d,nyb+(nguard-nyb/2)*k2d
           do i=1+nguard,nxb+nguard-nxb/2+1
             do ivar=1,nbndvar
               facevarx(ivar,i+ioff,j+joff,k+koff,jblock) = & 
     &                          send(ivar,i,j,k)
             enddo
           enddo
         enddo
       enddo

! y face next
       if(ndim.ge.2) then
       recv(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd) =  & 
     &  facevary(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd,isg)
       call amr_restrict_fc_fun(recv,temp,2)

       do k=1+nguard*k3d,nzb+nguard*k3d,2
         kk = (k-nguard*k3d)/2+1+nguard*k3d
         do j=1+nguard*k2d,nyb+(nguard+1)*k2d,2
           jj = (j-nguard)/2+1+nguard
           do i=1+nguard,nxb+nguard,2
             ii = (i-nguard)/2+1+nguard
             do ivar=1,nbndvar
               send(ivar,ii,jj,kk) = temp(ivar,i,j,k)
             enddo
           enddo
         enddo
       enddo

       do k=1+nguard*k3d,nzb+(nguard-nzb/2)*k3d
         do j=1+nguard*k2d,nyb+(nguard-nyb/2+1)*k2d
           do i=1+nguard,nxb+nguard-nxb/2
             do ivar=1,nbndvar
               facevary(ivar,i+ioff,j+joff,k+koff,jblock) = & 
     &                          send(ivar,i,j,k)
             enddo
           enddo
         enddo
       enddo
       endif

! z face last
       if(ndim.eq.3) then
       recv(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d) =  & 
     &  facevarz(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d, & 
     &  isg)
       call amr_restrict_fc_fun(recv,temp,3)

       do k=1+nguard*k3d,nzb+(nguard+1)*k3d,2
         kk = (k-nguard*k3d)/2+1+nguard*k3d
         do j=1+nguard*k2d,nyb+nguard*k2d,2
           jj = (j-nguard*k2d)/2+1+nguard*k2d
           do i=1+nguard,nxb+nguard,2
             ii = (i-nguard)/2+1+nguard
             do ivar=1,nbndvar
               send(ivar,ii,jj,kk) = temp(ivar,i,j,k)
             enddo
           enddo
         enddo
       enddo

       do k=1+nguard*k3d,nzb+(nguard-nzb/2+1)*k3d
         do j=1+nguard*k2d,nyb+(nguard-nyb/2)*k2d
           do i=1+nguard,nxb+nguard-nxb/2
             do ivar=1,nbndvar
               facevarz(ivar,i+ioff,j+joff,k+koff,jblock) = & 
     &                         send(ivar,i,j,k)
             enddo
           enddo
         enddo
       enddo
       endif

      endif

      endif
      enddo

! Now parents of leaf nodes which are on remote processors get data
! from their children and then perform restriction on it.

      do isg = 1,lnblocks

! Is this a parent block of at least one leaf node?
      if(nodetype(isg).eq.2.and.rflag(isg)) then

! If yes then cycle through its children.
      do ich=1,nchild

       jchild = ich

! Does the child reside on another processor?
       if(child(2,ich,isg).ne.mype) then

! Is this child a leaf block? If it is then fetch its data.
       remote_pe = child(2,ich,isg)
       remote_block  = child(1,ich,isg)
!       call shmem_integer_get(cnodetype,nodetype(remote_block),
!     .       1,remote_pe)
       cempty=1
!       call shmem_integer_get(cempty,empty(remote_block),
!     .       1,remote_pe)
       if(cnodetype.eq.1.and.cempty.eq.0) then

! compute the offset in the parent block appropriate for this child
       ioff = mod(jchild-1,2)*nxb/2
       joff = mod((jchild-1)/2,2)*nyb/2
       koff = mod((jchild-1)/4,2)*nzb/2



! x face first
!       call shmem_real_get
!     . ( recv(1,il_bnd,jl_bnd,
!     .       kl_bnd),
!     .       facevarx(1,1,1,1,remote_block),
!     .       len_blockfx*nfacevar,remote_pe)

! Compute restricted data from the data in the buffer
       call amr_restrict_fc_fun(recv,temp,1)

       do k=1+nguard*k3d,nzb+nguard*k3d,2
         kk = (k-nguard*k3d)/2+1+nguard*k3d
         do j=1+nguard*k2d,nyb+nguard*k2d,2
           jj = (j-nguard*k2d)/2+1+nguard*k2d
           do i=1+nguard,nxb+nguard+1,2
             ii = (i-nguard)/2+1+nguard
             do ivar=1,nbndvar
               send(ivar,ii,jj,kk) = temp(ivar,i,j,k)
             enddo
           enddo
         enddo
       enddo

! update the parent block
       do k=1+nguard*k3d,nzb+(nguard-nzb/2)*k3d
         do j=1+nguard*k2d,nyb+(nguard-nyb/2)*k2d
           do i=1+nguard,nxb+nguard-nxb/2+1
             do ivar=1,nbndvar
               facevarx(ivar,i+ioff,j+joff,k+koff,isg)= & 
     &                          send(ivar,i,j,k)
             enddo
           enddo
         enddo
       enddo



! y face next
       if(ndim.ge.2) then
!       call shmem_real_get
!     . ( recv(1,il_bnd,jl_bnd,
!     .       kl_bnd),
!     .       facevary(1,1,1,1,remote_block),
!     .       len_blockfy*nfacevar,remote_pe)

! Compute restricted data from the data in the buffer
       call amr_restrict_fc_fun(recv,temp,2)

       do k=1+nguard*k3d,nzb+nguard*k3d,2
         kk = (k-nguard*k3d)/2+1+nguard*k3d
         do j=1+nguard*k2d,nyb+(nguard+1)*k2d,2
           jj = (j-nguard*k2d)/2+1+nguard*k2d
           do i=1+nguard,nxb+nguard,2
             ii = (i-nguard)/2+1+nguard
             do ivar=1,nbndvar
               send(ivar,ii,jj,kk) = temp(ivar,i,j,k)
             enddo
           enddo
         enddo
       enddo

! update the parent block
       do k=1+nguard*k3d,nzb+(nguard-nzb/2)*k3d
         do j=1+nguard*k2d,nyb+(nguard-nyb/2+1)*k2d
           do i=1+nguard,nxb+nguard-nxb/2
             do ivar=1,nbndvar
               facevary(ivar,i+ioff,j+joff,k+koff,isg)= & 
     &                          send(ivar,i,j,k)
             enddo
           enddo
         enddo
       enddo
       endif

! z face last
       if(ndim.eq.3) then
!       call shmem_real_get
!     . ( recv(1,il_bnd,jl_bnd,
!     .       kl_bnd),
!     .       facevarz(1,1,1,1,remote_block),
!     .       len_blockfz*nfacevar,remote_pe)

! Compute restricted data from the data in the buffer
       call amr_restrict_fc_fun(recv,temp,3)

       do k=1+nguard*k3d,nzb+(nguard+1)*k3d,2
         kk = (k-nguard*k3d)/2+1+nguard*k3d
         do j=1+nguard*k2d,nyb+nguard*k2d,2
           jj = (j-nguard*k2d)/2+1+nguard*k2d
           do i=1+nguard,nxb+nguard,2
             ii = (i-nguard)/2+1+nguard
             do ivar=1,nbndvar
               send(ivar,ii,jj,kk) = temp(ivar,i,j,k)
             enddo
           enddo
         enddo
       enddo

! update the parent block
       do k=1+nguard*k3d,nzb+(nguard-nzb/2+1)*k3d
         do j=1+nguard*k2d,nyb+(nguard-nyb/2)*k2d
           do i=1+nguard,nxb+nguard-nxb/2
             do ivar=1,nbndvar
               facevarz(ivar,i+ioff,j+joff,k+koff,isg)= & 
     &                          send(ivar,i,j,k)
             enddo
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
!      call shmem_barrier_all()

      return
      end
