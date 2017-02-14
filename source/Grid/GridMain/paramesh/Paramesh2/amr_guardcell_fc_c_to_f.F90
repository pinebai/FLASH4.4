      subroutine amr_guardcell_fc_c_to_f(mype)


! $RCSfile: amr_guardcell_fc_c_to_f.F90,v $
! $Revision: 1.1 $
! $Date: 2004/05/20 01:03:44 $


!------------------------------------------------------------------------
! Rewritten to accomodate even block sizes but only 2nd order schemes
!
! This routine manages the transfer of face centered guard cell data 
! to all  leaf blocks from any neighbors which they have at a coarser 
! resolution.
!
! Written :     Peter MacNeice          July 1997
!------------------------------------------------------------------------
!
! Arguments:
!       mype            local processor number
!
!------------------------------------

use physicaldata
      use tree
      use workspace
      implicit none



      integer mype

      integer isg,jf,jblock,jchild,ich
      integer ioff,joff,koff

!------------------------------------
! local arrays
      real recv(nbndvar,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d, & 
     &       kl_bnd:ku_bnd+k3d)
      integer ichild(2,mchild)

      integer remote_pe,remote_block
      save recv

!------------------------------------

!      call shmem_barrier_all()

! cycle through the grid blocks on this processor
      if(lnblocks.gt.0) then
      do isg = 1,lnblocks

! Is this a leaf block and not at the original refinement level ?
      if( (nodetype(isg).eq.1) .and. & 
     &       lrefine(isg).gt.1) then


! Cycle over the blocks faces
       do jf = 1,nfaces

! Is the neighbor to this face at coarser resolution?
       if(neigh(1,jf,isg).eq.-1) then

! If yes then copy the appropriate layers from its parents data 
! into a temporary buffer array.


! identify which child this leaf block represents

! Does the parent reside on this processor?
       if(parent(2,isg).eq.mype) then
         jblock = parent(1,isg)
         ichild(:,:)=child(:,:,jblock)
       else
         remote_pe = parent(2,isg)
         remote_block  = parent(1,isg)
!         call shmem_integer_get(ichild,child(1,1,remote_block),
!     .       2*mchild,remote_pe)
       endif

       jchild = 0
       do ich=1,nchild
         if( (ichild(1,ich).eq.isg).and. & 
     &       (ichild(2,ich).eq.mype) ) jchild=ich
       enddo

! compute the offset in the parent block appropriate for this child
       ioff = mod(jchild-1,2)*nxb/2
       joff = mod((jchild-1)/2,2)*nyb/2
       koff = mod((jchild-1)/4,2)*nzb/2



! work on facevarx first ------------

! Does the parent reside on this processor?
       if(parent(2,isg).eq.mype) then

! copy required layers from parents data - note if we use the standard
! numerical recipes prolongation operator below, then we must copy 3 layers
! here, ie the guardcell layers and 2 immediate interior layers.

         recv(:,:,jl_bnd:ju_bnd,kl_bnd:ku_bnd)= & 
     &       facevarx(:,:,jl_bnd:ju_bnd,kl_bnd:ku_bnd,jblock)
       else
!         call shmem_real_get
!     . ( recv(1,1,
!     .       jl_bnd,kl_bnd),
!     .       facevarx(1,1,jl_bnd,
!     .       kl_bnd,remote_block),
!     .       len_blockfx*nbndvar,remote_pe)
       endif

! interpolate(prolongate) data from the parent to the child
       call amr_prolong_face_fun(recv,isg,ioff,joff,koff,jf,mype,1)



! now work on facevary  ------------
      if(ndim.ge.2) then

! Does the parent reside on this processor?
        if(parent(2,isg).eq.mype) then

! copy required layers from parents data - note if we use the standard
! numerical recipes prolongation operator below, then we must copy 3 layers
! here, ie the guardcell layers and 2 immediate interior layers.

           recv(:,il_bnd:iu_bnd,:,kl_bnd:ku_bnd)= & 
     &       facevary(:,il_bnd:iu_bnd,:, & 
     &       kl_bnd:ku_bnd,jblock)
        else
!           call shmem_real_get
!     . ( recv(1,
!     .       il_bnd,1,kl_bnd),
!     .       facevary(1,il_bnd,1,
!     .       kl_bnd,remote_block),
!     .       len_blockfy*nbndvar,remote_pe)
        endif

! interpolate(prolongate) data from the parent to the child
        call amr_prolong_face_fun(recv,isg,ioff,joff,koff,jf,mype,2)
      endif


! finally work on facevarz  ------------
      if(ndim.eq.3) then


! Does the parent reside on this processor?
        if(parent(2,isg).eq.mype) then

! copy required layers from parents data - note if we use the standard
! numerical recipes prolongation operator below, then we must copy 3 layers
! here, ie the guardcell layers and 2 immediate interior layers.

         recv(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,:)= & 
     &       facevarz(:,il_bnd:iu_bnd, & 
     &       jl_bnd:ju_bnd,:,jblock)
       else
!         call shmem_real_get
!     . ( recv(1,
!     .       il_bnd,jl_bnd,1),
!     .       facevarz(1,il_bnd,
!     .       jl_bnd,1,remote_block),
!     .       len_blockfz*nbndvar,remote_pe)

       endif

! interpolate(prolongate) data from the parent to the child
       call amr_prolong_face_fun(recv,isg,ioff,joff,koff,jf,mype,3)

      endif
! --------------------------------------


      endif

      enddo

      endif
      enddo
      endif

!      call shmem_barrier_all()


      return
      end
