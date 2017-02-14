      subroutine amr_face_cp_loc(idest,isource,iface,nlayers)


! $RCSfile: amr_face_cp_loc.F90,v $
! $Revision: 1.1 $
! $Date: 2004/05/20 01:03:44 $


!------------------------------------------------------------------------
!
! This routine copies guard cell information for data defined at
! cell centers to face iface in block
! idest, from the appropriate face of the neighboring block, assuming
! that both blocks are on the same processor.
! It can be easily edited to alter the data pattern required for schemes
! of different order.
!
! Arguments:
!      idest id of the destination block which requires guard 
!       cell info
!      isource id of the block from which guard cell info is to
!       be obtained
!      iface the id of the face of the block which requires
!       guard cell info
!       nlayers         the number of guard cell layers at each boundary
!
!
! Written :     Peter MacNeice          July 1997
!------------------------------------------------------------------------

use physicaldata
      use workspace
      implicit none

      integer idest,isource,iface,nlayers
      integer i_dest,i_source,j_dest,j_source,k_dest,k_source
      integer i,j,k



      if(iface.eq.2) then

       i_dest   = nxb+1+nguard + iface_off
       i_source = 1+nguard+gc_off_x + iface_off

       do i=0,nlayers-1
            facevarx(:,i_dest+i+1,:,:,idest) = & 
     &       facevarx(:,i_source+i+1,:,:,isource)
       enddo

       if(ndim.ge.2) then
       do i=0,nlayers-1
            facevary(:,i_dest+i,:,:,idest) = & 
     &       facevary(:,i_source+i,:,:,isource)
       enddo
       endif

       if(ndim.eq.3) then
       do i=0,nlayers-1
            facevarz(:,i_dest+i,:,:,idest) = & 
     &       facevarz(:,i_source+i,:,:,isource)
       enddo
       endif

      elseif(iface.eq.1) then

       i_dest   = 0+nguard + iface_off
       i_source = nxb+nguard-gc_off_x + iface_off

       do i=0,nlayers-1
            facevarx(:,i_dest-i,:,:,idest) = & 
     &       facevarx(:,i_source-i,:,:,isource)
       enddo

       if(ndim.ge.2) then
       do i=0,nlayers-1
            facevary(:,i_dest-i,:,:,idest) = & 
     &       facevary(:,i_source-i,:,:,isource)
       enddo
       endif

       if(ndim.eq.3) then
       do i=0,nlayers-1
            facevarz(:,i_dest-i,:,:,idest) = & 
     &       facevarz(:,i_source-i,:,:,isource)
       enddo
       endif

      elseif(iface.eq.4) then

       j_dest   = nyb+1+nguard + iface_off
       j_source = 1+nguard+gc_off_y + iface_off

       do j=0,nlayers-1
            facevarx(:,:,j_dest+j,:,idest) = & 
     &       facevarx(:,:,j_source+j,:,isource)
       enddo

       do j=0,nlayers-1
            facevary(:,:,j_dest+j+1,:,idest) = & 
     &       facevary(:,:,j_source+j+1,:,isource)
       enddo

       if(ndim.eq.3) then
       do j=0,nlayers-1
            facevarz(:,:,j_dest+j,:,idest) = & 
     &       facevarz(:,:,j_source+j,:,isource)
       enddo
       endif

      elseif(iface.eq.3) then

       j_dest   = 0+nguard + iface_off
       j_source = nyb+nguard-gc_off_y + iface_off

       do j=0,nlayers-1
            facevarx(:,:,j_dest-j,:,idest) = & 
     &       facevarx(:,:,j_source-j,:,isource)
       enddo

       do j=0,nlayers-1
            facevary(:,:,j_dest-j,:,idest) = & 
     &       facevary(:,:,j_source-j,:,isource)
       enddo

       if(ndim.eq.3) then
       do j=0,nlayers-1
            facevarz(:,:,j_dest-j,:,idest) = & 
     &       facevarz(:,:,j_source-j,:,isource)
       enddo
       endif

      elseif(iface.eq.6) then

       k_dest   = nzb+1+nguard + iface_off
       k_source = 1+nguard+gc_off_z + iface_off

       do k=0,nlayers-1
            facevarx(:,:,:,k_dest+k,idest) = & 
     &       facevarx(:,:,:,k_source+k,isource)
       enddo

       do k=0,nlayers-1
            facevary(:,:,:,k_dest+k,idest) = & 
     &       facevary(:,:,:,k_source+k,isource)
       enddo

       do k=0,nlayers-1
            facevarz(:,:,:,k_dest+k+1,idest) = & 
     &       facevarz(:,:,:,k_source+k+1,isource)
       enddo

      elseif(iface.eq.5) then

       k_dest   = 0+nguard + iface_off
       k_source = nzb+nguard-gc_off_z + iface_off

       do k=0,nlayers-1
            facevarx(:,:,:,k_dest-k,idest) = & 
     &       facevarx(:,:,:,k_source-k,isource)
       enddo

       do k=0,nlayers-1
            facevary(:,:,:,k_dest-k,idest) = & 
     &       facevary(:,:,:,k_source-k,isource)
       enddo

       do k=0,nlayers-1
            facevarz(:,:,:,k_dest-k,idest) = & 
     &       facevarz(:,:,:,k_source-k,isource)
       enddo

      endif

      return
      end
