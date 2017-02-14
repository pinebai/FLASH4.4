      subroutine amr_face_cp_remote(mype,remote_pe,remote_block,idest, & 
     &       iface,nlayers)


! $RCSfile: amr_face_cp_remote.F90,v $
! $Revision: 1.1 $
! $Date: 2004/05/20 01:03:44 $


!------------------------------------------------------------------------
!
! This routine copies guard cell information for cell face centered
! data to face iface in block
! idest, from the appropriate face of the neighboring block, assuming
! that the neighboring block is on a different processor.
! It can be easily edited to alter the data pattern required for schemes
! of different order.
!
! Arguments:
!      mype local processor
!      remote_pe remote processor
!      remote_block local block id of the block to be copied from
!       the remote processor
!      idest id of the destination block which requires guard 
!       cell info
!      iface the id of the face of the block which requires
!       guard cell info
!       nlayers         the number of guard cell layers at each boundary
!
!
!
! Written :     Peter MacNeice          July 1997
!------------------------------------------------------------------------

use physicaldata
      use tree
      use workspace
      implicit none


!-------------------------

      integer mype
      integer remote_pe,remote_block
      integer idest,iface,nlayers

! local arrays
      real recvx(nbndvar,il_bnd:iu_bnd+1,jl_bnd:ju_bnd, & 
     &       kl_bnd:ku_bnd)
      real recvy(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d, & 
     &       kl_bnd:ku_bnd)
      real recvz(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd, & 
     &       kl_bnd:ku_bnd+k3d)

      integer i,j,k
      integer i_dest,i_source,j_dest,j_source,k_dest,k_source

      save recvx,recvy,recvz

!-------------------------


! Copy complete remote block into a buffer block called recv.
!      call shmem_real_get
!     . (recvx,facevarx(1,1,1,1,remote_block),
!     .       len_blockfx*nbndvar,remote_pe) 
!      call shmem_real_get
!     . (recvy,facevary(1,1,1,1,remote_block),
!     .       len_blockfy*nbndvar,remote_pe) 
!      call shmem_real_get
!     . (recvz,facevarz(1,1,1,1,remote_block),
!     .       len_blockfz*nbndvar,remote_pe) 


! Select the correct face data out of the recv block.
      if(iface.eq.2) then

       i_dest   = nxb+1+nguard + iface_off
       i_source = 1+nguard+gc_off_x + iface_off

       do i=0,nlayers-1
         facevarx(:,i_dest+i+1,:,:,idest)= recvx(:,i_source+i+1,:,:)
       enddo

       if(ndim.ge.2) then
            do i=0,nlayers-1
              facevary(:,i_dest+i,:,:,idest)= recvy(:,i_source+i,:,:)
            enddo
       endif

       if(ndim.eq.3) then
            do i=0,nlayers-1
              facevarz(:,i_dest+i,:,:,idest)= recvz(:,i_source+i,:,:)
            enddo
       endif

      elseif(iface.eq.1) then

       i_dest   = 0+nguard + iface_off
       i_source = nxb+nguard-gc_off_x + iface_off

       do i=0,nlayers-1
         facevarx(:,i_dest-i,:,:,idest)= recvx(:,i_source-i,:,:)
       enddo

       if(ndim.ge.2) then
         do i=0,nlayers-1
           facevary(:,i_dest-i,:,:,idest)= recvy(:,i_source-i,:,:)
         enddo
       endif

       if(ndim.eq.3) then
         do i=0,nlayers-1
           facevarz(:,i_dest-i,:,:,idest)= recvz(:,i_source-i,:,:)
         enddo
       endif

      elseif(iface.eq.4) then

       j_dest   = nyb+1+nguard + iface_off
       j_source = 1+nguard+gc_off_y + iface_off

                do j=0,nlayers-1
            facevarx(:,:,j_dest+j,:,idest)= & 
     &       recvx(:,:,j_source+j,:)
       enddo

                do j=0,nlayers-1
            facevary(:,:,j_dest+j+1,:,idest)= & 
     &       recvy(:,:,j_source+j+1,:)
       enddo

       if(ndim.eq.3) then
                do j=0,nlayers-1
            facevarz(:,:,j_dest+j,:,idest)= & 
     &       recvz(:,:,j_source+j,:)
       enddo
       endif

      elseif(iface.eq.3) then

       j_dest   = 0+nguard + iface_off
       j_source = nyb+nguard-gc_off_y + iface_off

                do j=0,nlayers-1
            facevarx(:,:,j_dest-j,:,idest)= & 
     &       recvx(:,:,j_source-j,:)
       enddo

                do j=0,nlayers-1
            facevary(:,:,j_dest-j,:,idest)= & 
     &       recvy(:,:,j_source-j,:)
       enddo

       if(ndim.eq.3) then
                do j=0,nlayers-1
            facevarz(:,:,j_dest-j,:,idest)= & 
     &       recvz(:,:,j_source-j,:)
       enddo
       endif

      elseif(iface.eq.6) then

       k_dest   = nzb+1+nguard + iface_off
       k_source = 1+nguard+gc_off_z + iface_off

                do k=0,nlayers-1
            facevarx(:,:,:,k_dest+k,idest)= & 
     &       recvx(:,:,:,k_source+k)
       enddo

                do k=0,nlayers-1
            facevary(:,:,:,k_dest+k,idest)= & 
     &       recvy(:,:,:,k_source+k)
       enddo

                do k=0,nlayers-1
            facevarz(:,:,:,k_dest+k+1,idest)= & 
     &       recvz(:,:,:,k_source+k+1)
       enddo

      elseif(iface.eq.5) then

       k_dest   = 0+nguard + iface_off
       k_source = nzb+nguard-gc_off_z + iface_off

                do k=0,nlayers-1
            facevarx(:,:,:,k_dest-k,idest)= & 
     &       recvx(:,:,:,k_source-k)
       enddo

                do k=0,nlayers-1
            facevary(:,:,:,k_dest-k,idest)= & 
     &       recvy(:,:,:,k_source-k)
       enddo

                do k=0,nlayers-1
            facevarz(:,:,:,k_dest-k,idest)= & 
     &       recvz(:,:,:,k_source-k)
       enddo

      endif

      return
      end
