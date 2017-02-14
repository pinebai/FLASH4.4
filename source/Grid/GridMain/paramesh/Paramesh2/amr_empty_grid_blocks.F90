      subroutine amr_empty_grid_blocks(mype,iopt)


! $RCSfile: amr_empty_grid_blocks.F90,v $
! $Revision: 1.1 $
! $Date: 2004/05/20 01:03:44 $


!------------------------------------------------------------------------
!
! This routine sets up any intentionally empty interior grid blocks
! with the data required to ensure that the surrounding active
! blocks get the correct data during guard cell filling. It does this
! by applying a boundary condition routine on the appropriate faces
! of the surrounding blocks, and then copying this data into the
! equivalent interior cells of the empty block.
!
! Note, empty blocks are assigned empty=1.
!
! Written :     Peter MacNeice          August 1997
!------------------------------------------------------------------------

use physicaldata
      use tree
      use workspace
      implicit none



!------------------------------------
! local variables

      real recv(nvar,1:iu_bnd,1:ju_bnd,1:ku_bnd)
      real recvfx(nbndvar,1:iu_bnd+1,1:ju_bnd,1:ku_bnd)
      real recvfy(nbndvar,1:iu_bnd,1:ju_bnd+k2d,1:ku_bnd)
      real recvfz(nbndvar,1:iu_bnd,1:ju_bnd,1:ku_bnd+k3d)

      integer remote_pe,remote_block
      integer parent_pe,parent_block
      integer cempty,cnodetype,cneigh(2),ichild(2,mchild)
      save    cempty,cnodetype,cneigh,ichild

      save    recv,recvfx,recvfy,recvfz

!------------------------------------

      parameter(lengthc_x = nvar*nguard*(nyb+2*nguard*k2d) & 
     &       *(nzb+2*nguard*k3d))
      parameter(lengthc_y = nvar*(nxb+2*nguard)*nguard & 
     &       *(nzb+2*nguard*k3d))
      parameter(lengthc_z = nvar*(nxb+2*nguard) & 
     &       *(nyb+2*nguard*k2d)*nguard)
      parameter(lengthw_x = nguard_work*(nyb+2*nguard_work*k2d) & 
     &       *(nzb+2*nguard_work*k3d))
      parameter(lengthw_y = (nxb+2*nguard_work)*nguard_work & 
     &       *(nzb+2*nguard_work*k3d))
      parameter(lengthw_z = (nxb+2*nguard_work) & 
     &       *(nyb+2*nguard_work*k2d)*nguard_work)
      parameter(lengthf_x = nbndvar*(nguard+1)*(nyb+2*nguard*k2d) & 
     &       *(nzb+2*nguard*k3d))
      parameter(lengthf_y = nbndvar*(nxb+2*nguard)*(nguard+k2d) & 
     &       *(nzb+2*nguard*k3d))
      parameter(lengthf_z = nbndvar*(nxb+2*nguard)* & 
     &       (nyb+2*nguard*k2d)*(nguard+k3d))


      integer ka,kb,kc,kd
      parameter(ka=1+nguard*k3d)
      parameter(kb=1+(2*nguard-1)*k3d)
      parameter(kc=1+k3d*(nzb+nguard))
      parameter(kd=nzb+2*nguard*k3d)

!----------------------------------------------------------------
! Step 1.   Set boundary layers on faces of leaf blocks which face
!           empty leaf blocks or the parents of empty leaf blocks.


! cycle through the grid blocks on this processor
      if(lnblocks.gt.0) then
      do isg = 1,lnblocks

! Is this a leaf block which is not empty?
      if(nodetype(isg).eq.1.and.empty(isg).eq.0) then


! Cycle over the blocks faces
       do jf = 1,nfaces

            remote_pe = neigh(2,jf,isg)
            remote_block  = neigh(1,jf,isg)

! Is the neighbor to this face an empty leaf block or the parent of
! an empty leaf block?
! If yes then set the boundary guard cells on this block.
            cempty = 0
!            if(remote_block.gt.0) call shmem_integer_get(cempty,
!     .                             empty(remote_block),1,remote_pe)
            if(cempty.eq.1.or.cempty.eq.2) then
                   call amr_bc_block(jf,remote_block,iopt,isg,mype)
            endif

       enddo


      endif
      enddo
      endif

!      call shmem_barrier_all()

!----------------------------------------------------------------
! Step 2.   Set boundary layers on faces of the parents of leaf 
!           blocks if those parents face empty leaf blocks.



! cycle through the grid blocks on this processor
        if(lnblocks.gt.0) then
        do isg = 1,lnblocks

! Is this a parent of a leaf block which is not empty?
        if(nodetype(isg).eq.2.and.empty(isg).eq.0) then


! Cycle over the blocks faces
        do jf = 1,nfaces

            remote_pe = neigh(2,jf,isg)
            remote_block  = neigh(1,jf,isg)

! Is the neighbor to this face an empty leaf block?
! If yes then set the boundary guard cells on this block.
            cempty = 0
!            if(remote_block.gt.0) call shmem_integer_get(cempty,
!     .                          empty(remote_block),1,remote_pe)
            if(cempty.eq.1.or.cempty.eq.2) then
                    call amr_bc_block(jf,remote_block,iopt,isg,mype)
            endif

        enddo

        endif
        enddo
        endif

!      call shmem_barrier_all()


!----------------------------------------------------------------
! Step 3.    Now do SRL filling of the interior layers next to the 
!            boundaries of empty blocks.


! cycle through the grid blocks on this processor
        if(lnblocks.gt.0) then
        do isg = 1,lnblocks

! Is this an empty leaf block or a parent of an empty leaf block ?
        if(empty(isg).eq.1.or.empty(isg).eq.2) then


! Cycle over the blocks faces
        do jf = 1,nfaces

            remote_pe = neigh(2,jf,isg)
            remote_block  = neigh(1,jf,isg)


! If a non-empty neighbor exists at this refinement level copy its data 
! into temporary local storage.
            if(remote_block.gt.0) then

                cempty = 1
!                call shmem_integer_get(cempty,empty(remote_block),
!     .                                                 1,remote_pe)
                if(cempty.eq.0) then

                  if(iopt.eq.1) then
!                      call shmem_real_get
!     . (recv,
!     .                 unk(1,1,1,1,remote_block),len_block,remote_pe)

                      if(nfacevar.gt.0) then
!                          call shmem_real_get
!     . (recvfx,
!     .                                facevarx(1,1,1,1,remote_block),
!     .                                len_blockfx*nbndvar,remote_pe)
!                          call shmem_real_get
!     . (recvfy,
!     .                                facevary(1,1,1,1,remote_block),
!     .                                len_blockfy*nbndvar,remote_pe)
!                          if(ndim.eq.3) call shmem_real_get
!     . (recvfz,
!     .                                facevarz(1,1,1,1,remote_block),
!     .                                len_blockfz*nbndvar,remote_pe)
                      endif


                      if(jf.eq.1) then

                       unk(:,1+nguard:2*nguard,:,:,isg)= & 
     &                       recv(:,1+nxb+nguard:nxb+2*nguard,:,:)
                       if(nfacevar.gt.0) then
                         facevarx(:,1+nguard:2*nguard+1,:,:,isg)= & 
     &                     recvfx(:,1+nxb+nguard:1+nxb+2*nguard,:,:)
                         facevary(:,1+nguard:2*nguard,:,:,isg)= & 
     &                     recvfy(:,1+nxb+nguard:nxb+2*nguard,:,:)
                         if(ndim.eq.3)  & 
     &                     facevarz(:,1+nguard:2*nguard,:,:,isg)= & 
     &                     recvfz(:,1+nxb+nguard:nxb+2*nguard,:,:)
                       endif

                      elseif(jf.eq.2) then

                       unk(:,1+nxb:nxb+nguard,:,:,isg)= & 
     &                       recv(:,1:nguard,:,:)
                       if(nfacevar.gt.0) then
                         facevarx(:,1+nxb:1+nxb+nguard,:,:,isg) & 
     &                     =recvfx(:,1:1+nguard,:,:)
                         facevary(:,1+nxb:nxb+nguard,:,:,isg) & 
     &                     =recvfy(:,1:nguard,:,:)
                       if(ndim.eq.3)  & 
     &                   facevarz(:,1+nxb:nxb+nguard,:,:,isg) & 
     &                   =recvfz(:,1:nguard,:,:)
                       endif

                      elseif(jf.eq.3) then

                       unk(:,:,1+nguard:2*nguard,:,isg)= & 
     &                       recv(:,:,1+nyb+nguard:nyb+2*nguard,:)
                       if(nfacevar.gt.0) then
                         facevary(:,:,1+nguard:1+2*nguard,:,isg) & 
     &                     =recvfy(:,:,1+nyb+nguard:1+nyb+2*nguard,:)
                         facevarx(:,:,1+nguard:2*nguard,:,isg) & 
     &                     =recvfx(:,:,1+nyb+nguard:nyb+2*nguard,:)
                         if(ndim.eq.3)  & 
     &                     facevarz(:,:,1+nguard:2*nguard,:,isg) & 
     &                     =recvfz(:,:,1+nyb+nguard:nyb+2*nguard,:)
                       endif


                      elseif(jf.eq.4) then

                       unk(:,:,1+nyb:nyb+nguard,:,isg)= & 
     &                       recv(:,:,1:nguard,:)
                       if(nfacevar.gt.0) then
                         facevary(:,:,1+nyb:1+nyb+nguard,:,isg) & 
     &                     =recvfy(:,:,1:1+nguard,:)
                         facevarx(:,:,1+nyb:nyb+nguard,:,isg) & 
     &                     =recvfx(:,:,1:nguard,:)
                         if(ndim.eq.3)  & 
     &                     facevarz(:,:,1+nyb:nyb+nguard,:,isg) & 
     &                     =recvfz(:,:,1:nguard,:)
                       endif


                      elseif(jf.eq.5) then

                       unk(:,:,:,ka:kb,isg) = recv(:,:,:,kc:kd)
                       if(nfacevar.gt.0) then
                         facevarz(:,:,:,ka:kb+k3d,isg) =  & 
     &                     recvfz(:,:,:,kc:kd+k3d)
                         facevarx(:,:,:,ka:kb,isg) =  & 
     &                     recvfx(:,:,:,kc:kd)
                         facevary(:,:,:,ka:kb,isg) =  & 
     &                     recvfy(:,:,:,kc:kd)
                       endif


                      elseif(jf.eq.6) then

                       unk(:,:,:,k3d+nzb:nzb+k3d*nguard,isg)= & 
     &                       recv(:,:,:,1:1+k3d*(nguard-1))
                       if(nfacevar.gt.0) then
                          facevarz(:,:,:,k3d+nzb:kc,isg) & 
     &                      =recvfz(:,:,:,1:ka)
                          facevarx(:,:,:,k3d+nzb:nzb+k3d*nguard,isg) & 
     &                      =recvfx(:,:,:,1:1+k3d*(nguard-1))
                          facevary(:,:,:,k3d+nzb:nzb+k3d*nguard,isg) & 
     &                      =recvfy(:,:,:,1:1+k3d*(nguard-1))
                       endif

                      endif




       elseif(iopt.eq.2) then

!       call shmem_real_get
!     . (recv1,work(1,1,1,remote_block),
!     .       len_wblock,remote_pe)

       if(jf.eq.1) then
            work(1+nguard_work:2*nguard_work,:,:,isg)= & 
     &       recv1(1+nxb+nguard_work:nxb+2*nguard_work,:,:)
       elseif(jf.eq.2) then
            work(1+nxb:nxb+nguard_work,:,:,isg)= & 
     &       recv1(1:nguard_work,:,:)
       elseif(jf.eq.3) then
            work(:,1+nguard_work:2*nguard_work,:,isg)= & 
     &       recv1(:,1+nyb+nguard_work:nyb+2*nguard_work,:)
       elseif(jf.eq.4) then
                        work(:,1+nyb:nyb+nguard_work,:,isg)= & 
     &                  recv1(:,1:nguard_work,:)
       elseif(jf.eq.5) then
            work(:,:,1+nguard_work*k3d:1+k3d*(2*nguard_work-1),isg)= & 
     &      recv1(:,:,1+k3d*(nzb+nguard_work):nzb+2*nguard_work*k3d)

       elseif(jf.eq.6) then
                        work(:,:,k3d+nzb:nzb+k3d*nguard_work,isg)= & 
     &                  recv1(:,:,1:1+k3d*(nguard_work-1))
       endif


       endif
       endif                              ! test if neigh is not empty
       endif                              ! test if neigh exists

       enddo                              ! end of loop over faces

      endif
      enddo
      endif

!      call shmem_barrier_all()

!----------------------------------------------------------------


      return
      end
