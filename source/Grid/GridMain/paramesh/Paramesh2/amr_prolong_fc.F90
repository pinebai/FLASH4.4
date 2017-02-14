      subroutine amr_prolong_fc(mype)


! $RCSfile: amr_prolong_fc.F90,v $
! $Revision: 1.1 $
! $Date: 2004/05/20 01:03:44 $


!------------------------------------------------------------------------
!
! This routine interpolates data from a parent block to any
! newly created child blocks, for data defined at cell face centers.
!
! This routine calls user supplied routines called prolong_face_fun
! which contains the detailed definition of the prolongation operator
! which is to be applied to the data in facevarx(y)(z).
!
! Written :     Peter MacNeice          July 1997
!------------------------------------------------------------------------

use physicaldata
      use tree
      use workspace
      implicit none


!------------------------------------
! local arrays
      real recv(nbndvar,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d, & 
     &       kl_bnd:ku_bnd+k3d)
      integer ichild(2,mchild)

      integer remote_pe,remote_block
      integer jchild,jblock,jface
      integer mype,ich
      integer ioff,joff,koff

      integer isg

      save recv

!------------------------------------


      if(lnblocks.gt.0) then


! cycle through the grid blocks on this processor
      do isg = 1,lnblocks

! Is this a new child block ?
      if(newchild(isg)) then

! Does the parent reside on this processor?
      if(parent(2,isg).eq.mype) then

        jblock = parent(1,isg)

! identify which child this leaf block represents
        jchild = 0
        do ich=1,nchild
          if( (child(1,ich,parent(1,isg)).eq.isg) & 
     &       .and.(child(2,ich,parent(1,isg)).eq.mype) ) & 
     &       jchild=ich
        enddo

! compute the offset in the parent block appropriate for this child
        ioff = mod(jchild-1,2)*nxb/2
        joff = mod((jchild-1)/2,2)*nyb/2
        koff = mod((jchild-1)/4,2)*nzb/2

! select complete block for prolongation
        jface=0


! x-face first
! copy its parents data into a temporary buffer array.
        recv(:,:,jl_bnd:ju_bnd,kl_bnd:ku_bnd) = & 
     &       facevarx(:,:,jl_bnd:ju_bnd,kl_bnd:ku_bnd,jblock)

! interpolate data from the parent to the child
        call amr_prolong_face_fun(recv,isg,ioff,joff,koff,jface,mype,1)

! y-face
! copy its parents data into a temporary buffer array.
        if(ndim.ge.2) then
          recv(:,il_bnd:iu_bnd,:,kl_bnd:ku_bnd) = & 
     &       facevary(:,il_bnd:iu_bnd,:,kl_bnd:ku_bnd,jblock)

! interpolate data from the parent to the child
          call amr_prolong_face_fun(recv,isg,ioff,joff,koff,jface, & 
     &                          mype,2)
        endif

! z-face
! copy its parents data into a temporary buffer array.
        if(ndim.eq.3) then
          recv(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,:) = & 
     &       facevarz(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,:,jblock)

! interpolate data from the parent to the child
          call amr_prolong_face_fun(recv,isg,ioff,joff,koff,jface, & 
     &                          mype,3)
        endif


      endif

      endif
      enddo


! At this point every new child block whose parent is on this processor
! has received its data. It remains for new children whose parents are off
! processor to get their data.


      do isg = 1,lnblocks

! Is this a new child block ?
      if(newchild(isg)) then

! Does the parent reside on this processor?
      if(parent(2,isg).ne.mype) then

        remote_pe = parent(2,isg)
        remote_block  = parent(1,isg)
!        call shmem_integer_get(ichild,child(1,1,remote_block),
!     .       2*mchild,remote_pe)

! identify which child this leaf block represents
        jchild = 0
        do ich=1,nchild
          if( (ichild(1,ich).eq.isg).and. & 
     &       (ichild(2,ich).eq.mype) ) jchild=ich
        enddo

! compute the offset in the parent block appropriate for this child
        ioff = mod(jchild-1,2)*nxb/2
        joff = mod((jchild-1)/2,2)*nyb/2
        koff = mod((jchild-1)/4,2)*nzb/2

! interpolate data from the parent to the child
        jface=0


! x-face
!        call shmem_real_get
!     . (recv(1,1,jl_bnd,kl_bnd),
!     .       facevarx(1,1,jl_bnd,kl_bnd,remote_block),
!     .       len_blockfx*nbndvar,remote_pe)
        call amr_prolong_face_fun(recv,isg,ioff,joff,koff,jface,mype,1)

! y-face
        if(ndim.ge.2) then
!          call shmem_real_get
!     . (recv(1,il_bnd,1,kl_bnd),
!     .       facevary(1,il_bnd,1,kl_bnd,remote_block),
!     .       len_blockfy*nbndvar,remote_pe)
          call amr_prolong_face_fun(recv,isg,ioff,joff,koff,jface, & 
     &                          mype,2)
        endif

! z-face
        if(ndim.eq.3) then
!          call shmem_real_get
!     . (recv(1,il_bnd,jl_bnd,1),
!     .       facevarz(1,il_bnd,jl_bnd,1,remote_block),
!     .       len_blockfz*nbndvar,remote_pe)
          call amr_prolong_face_fun(recv,isg,ioff,joff,koff,jface, & 
     &                          mype,3)
        endif

        endif

        endif
        enddo


      endif

!      call shmem_barrier_all()

! Ensure new blocks inherit data on block face shared with an old
! existing neighbor, instead of filling from parent by interpolation.
      call amr_prolong_fc_divbconsist(mype)

!      call shmem_barrier_all()

      return
      end

