      subroutine amr_restrict_edge_data(mype)


! $RCSfile: amr_restrict_edge_data.F90,v $
! $Revision: 1.1 $
! $Date: 2004/05/20 01:03:44 $


!------------------------------------------------------------------------
!
! This routine does the data averaging on cell edges required when a 
! child block passes data back to its parent. The parent receives data 
! at the block boundary only.
!
! This routine provides a mechanism for passing data defined at block
! boundaries from leaf blocks back to their parents.
! The averaging rules used to combine interface values on the finer
! mesh to construct interface values on the coarser parent mesh are
! specified by the user who provides a function called amr_restrict_edge
! to do this.
!
! This routine is only relevant for schemes with even number of grid points.
!
!
! Written :     Peter MacNeice          July 1997
!------------------------------------------------------------------------

use physicaldata
      use tree
      implicit none



!------------------------------------
! local arrays

      integer remote_pe,remote_block
      integer cnodetype

      integer mype
      integer ich,isg,icoord
      integer jchild,jblock
      integer iface,jface,kface
      integer ioff,joff,koff
      integer i,j,k,ii,jj,kk
      
!------------------------------------


! cycle through the grid blocks on this processor
      if(lnblocks.gt.0) then
      do isg = 1,lnblocks

! Is this a parent block of at least one leaf node?
      if(nodetype(isg).eq.2) then

! If yes then cycle through its children.
       do ich=1,nchild

       jchild = ich
       remote_pe = child(2,ich,isg)
       remote_block  = child(1,ich,isg)

! Is this child a leaf block(nodetype=1)? 
! If it is then fetch its data.
!       call shmem_integer_get(cnodetype,nodetype(remote_block),
!     .       1,remote_pe)
       if(cnodetype.eq.1) then

       do icoord=1,ndim
!       if(icoord.eq.1) then
!       call shmem_real_get
!     . (recvarx1e,
!     .       bedge_facex_y(1,1,1,1,remote_block),
!     .       len_block_ex*nedges,remote_pe)
!       if((ndim.eq.3).or.(l2p5d.eq.1))
!     .       call shmem_real_get
!     . (recvarx2e,
!     .       bedge_facex_z(1,1,1,1,remote_block),
!     .       len_block_ex*nedges,remote_pe)
!       elseif(icoord.eq.2) then
!       if((ndim.eq.3).or.(l2p5d.eq.1))
!     .       call shmem_real_get
!     . (recvary1e,
!     .       bedge_facey_z(1,1,1,1,remote_block),
!     .       len_block_ey*nedges,remote_pe)
!       call shmem_real_get
!     . (recvary2e,
!     .       bedge_facey_x(1,1,1,1,remote_block),
!     .       len_block_ey*nedges,remote_pe)
!       elseif(icoord.eq.3) then
!       call shmem_real_get
!     . (recvarz1e,
!     .       bedge_facez_x(1,1,1,1,remote_block),
!     .       len_block_ez*nedges,remote_pe)
!       call shmem_real_get
!     . (recvarz2e,
!     .       bedge_facez_y(1,1,1,1,remote_block),
!     .       len_block_ez*nedges,remote_pe)
!       endif


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

       call amr_restrict_edge(icoord)
       i = iface+1
                        do k=1+nguard*k3d,nzb+(nguard+1)*k3d,2
       kk = (k-nguard*k3d)/2+nguard*k3d+1
                        do j=1+nguard,nyb+nguard,2
       jj = (j-nguard)/2+nguard+1
                        bedge_facex_y(:,i,jj+joff,kk+koff, & 
     &       isg) = recvarx1e(:,i,j,k)
                        enddo
                        enddo

       if((ndim.eq.3).or.(l2p5d.eq.1)) then
                        do k=1+nguard*k3d,nzb+nguard*k3d,2
       kk = (k-nguard*k3d)/2+nguard*k3d+1
                        do j=1+nguard,nyb+nguard+1,2
       jj = (j-nguard)/2+nguard+1
                        bedge_facex_z(:,i,jj+joff,kk+koff, & 
     &       isg) = recvarx2e(:,i,j,k)
                        enddo
                        enddo
       endif

                        elseif(icoord.eq.2) then

       call amr_restrict_edge(icoord)
       j = jface+1
                        do k=1+nguard*k3d,nzb+(nguard+1)*k3d,2
       kk = (k-nguard*k3d)/2+nguard*k3d+1
                        do i=1+nguard,nxb+nguard,2
       ii = (i-nguard)/2+nguard+1
                        bedge_facey_x(:,ii+ioff,j,kk+koff, & 
     &       isg) = recvary2e(:,i,j,k)

                        enddo
                        enddo

       if((ndim.eq.3).or.(l2p5d.eq.1)) then
       do k=1+nguard*k3d,nzb+nguard*k3d,2
       kk = (k-nguard*k3d)/2+nguard*k3d+1
                        do i=1+nguard,nxb+nguard+1,2
       ii = (i-nguard)/2+nguard+1
                        bedge_facey_z(:,ii+ioff,j,kk+koff, & 
     &       isg) = recvary1e(:,i,j,k)
                        enddo
                        enddo
       endif

                        elseif(icoord.eq.3) then

       call amr_restrict_edge(icoord)
       k = kface+1
                        do j=1+nguard,nyb+nguard+1,2
       jj = (j-nguard)/2+nguard+1
                        do i=1+nguard,nxb+nguard,2
       ii = (i-nguard)/2+nguard+1
                        bedge_facez_x(:,ii+ioff,jj+joff,k, & 
     &       isg) = recvarz1e(:,i,j,k)
                        enddo
                        enddo
                        do j=1+nguard,nyb+nguard,2
       jj = (j-nguard)/2+nguard+1
                        do i=1+nguard,nxb+nguard+1,2
       ii = (i-nguard)/2+nguard+1
                        bedge_facez_y(:,ii+ioff,jj+joff,k, & 
     &       isg) = recvarz2e(:,i,j,k)
                        enddo
                        enddo

                        endif

       enddo

       endif

       enddo

      endif

      enddo
      endif

!      call shmem_barrier_all()

      return
      end
