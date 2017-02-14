      subroutine amr_cp_remote(mype,remote_pe,remote_block,idest, & 
     &       iface,iopt,nlayers)


! $RCSfile: amr_cp_remote.F90,v $
! $Revision: 1.1 $
! $Date: 2004/05/20 01:03:44 $


!------------------------------------------------------------------------
!
! This routine copies guard cell information to face iface in block
! idest, from the appropriate face of the neighboring block, assuming
! that the neighboring block is on a different processor.
! It can be easily edited to alter the data pattern required for schemes
! of different order.
!
! Arguments:
!      mype             local processor
!      remote_pe        remote processor
!      remote_block     local block id of the block to be copied from
!                        the remote processor
!      idest            id of the destination block which requires guard 
!                        cell info
!      iface            the id of the face of the block which requires
!                        guard cell info
!      iopt             a switch to control which data source is to be used
!                        iopt=1 will use 'unk'
!                        iopt=2 will use 'work'
!      nlayers          the number of guard cell layers at each boundary
!
!
!
! Written :     Peter MacNeice          February 1997
!------------------------------------------------------------------------

use physicaldata
      use tree
      use workspace
      implicit none



#ifdef TIMINGS
#include "timer.fh"
#endif

!-------------------------

!      integer remote_pe,remote_block
      integer mype,remote_pe,remote_block,idest,iface,iopt,nlayers

      integer i_dest,i_source,j_dest,j_source,k_dest,k_source
      integer i,j,k,ivar

! local arrays
      real recv(nvar,il_bnd:iu_bnd,jl_bnd:ju_bnd, & 
     &       kl_bnd:ku_bnd)
      save recv

!-------------------------

#ifdef TIMINGS
      itimer1 = irtc()
#endif

! Copy complete remote block into a buffer block called recv.
      if(iopt.eq.1) then
!      call shmem_real_get
!     . (recv,unk(1,1,1,1,remote_block),len_block,
!     .       remote_pe) 
      elseif(iopt.eq.2) then
!      call shmem_real_get
!     . (recv1,work(1,1,1,remote_block),len_wblock,
!     .       remote_pe) 
      endif



!-------------------------
      if(iopt.eq.1) then
!-------------------------

! Select the correct face data out of the recv block.
      if(iface.eq.2) then

       i_dest   = nxb+1+nguard
       i_source = 1+nguard+gc_off_x
       
       do i=0,nlayers-1
       do ivar=1,nvar
            unk(ivar,i_dest+i,:,:,idest)= recv(ivar,i_source+i,:,:)
       enddo
       enddo

      elseif(iface.eq.1) then

       i_dest   = 0+nguard
       i_source = nxb+nguard-gc_off_x

       do i=0,nlayers-1
       do ivar=1,nvar
            unk(ivar,i_dest-i,:,:,idest)= recv(ivar,i_source-i,:,:)
       enddo
       enddo

      elseif(iface.eq.4) then

       j_dest   = nyb+1+nguard
       j_source = 1+nguard+gc_off_y

       do j=0,nlayers-1
       do ivar=1,nvar
            unk(ivar,:,j_dest+j,:,idest)= recv(ivar,:,j_source+j,:)
       enddo
       enddo

      elseif(iface.eq.3) then

       j_dest   = 0+nguard
       j_source = nyb+nguard-gc_off_y

       do j=0,nlayers-1
       do ivar=1,nvar
            unk(ivar,:,j_dest-j,:,idest)= recv(ivar,:,j_source-j,:)
       enddo
       enddo

      elseif(iface.eq.6) then

       k_dest   = nzb+1+nguard
       k_source = 1+nguard+gc_off_z

       do k=0,nlayers-1
       do ivar=1,nvar
            unk(ivar,:,:,k_dest+k,idest)= recv(ivar,:,:,k_source+k)
       enddo
       enddo

      elseif(iface.eq.5) then

       k_dest   = 0+nguard
       k_source = nzb+nguard-gc_off_z

       do k=0,nlayers-1
       do ivar=1,nvar
            unk(ivar,:,:,k_dest-k,idest)= recv(ivar,:,:,k_source-k)
       enddo
       enddo

      endif



!-------------------------
      elseif(iopt.eq.2) then
!-------------------------


! Select the correct face data out of the recv block.
      if(iface.eq.2) then

       i_dest   = nxb+1+nguard_work
       i_source = 1+nguard_work+gc_off_x

       do i=0,nlayers-1
            work(i_dest+i,:,:,idest,1)= recv1(i_source+i,:,:)
       enddo

      elseif(iface.eq.1) then

       i_dest   = 0+nguard_work
       i_source = nxb+nguard_work-gc_off_x

       do i=0,nlayers-1
            work(i_dest-i,:,:,idest,1)= recv1(i_source-i,:,:)
       enddo

      elseif(iface.eq.4) then

       j_dest   = nyb+1+nguard_work
       j_source = 1+nguard_work+gc_off_y

       do j=0,nlayers-1
            work(:,j_dest+j,:,idest,1)= recv1(:,j_source+j,:)
       enddo

      elseif(iface.eq.3) then

       j_dest   = 0+nguard_work
       j_source = nyb+nguard_work-gc_off_y

       do j=0,nlayers-1
            work(:,j_dest-j,:,idest,1)= recv1(:,j_source-j,:)
       enddo

      elseif(iface.eq.6) then

       k_dest   = nzb+1+nguard_work
       k_source = 1+nguard_work+gc_off_z

       do k=0,nlayers-1
            work(:,:,k_dest+k,idest,1)= recv1(:,:,k_source+k)
       enddo

      elseif(iface.eq.5) then

       k_dest   = 0+nguard_work
       k_source = nzb+nguard_work-gc_off_z

       do k=0,nlayers-1
            work(:,:,k_dest-k,idest,1)= recv1(:,:,k_source-k)
       enddo

      endif


!-------------------------
      endif
!-------------------------


#ifdef TIMINGS
      itimer2 = irtc()
      irtc_cprem = itimer2-itimer1+irtc_cprem
#endif

      return
      end
