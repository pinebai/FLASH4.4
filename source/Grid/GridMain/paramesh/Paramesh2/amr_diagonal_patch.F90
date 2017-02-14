      subroutine amr_diagonal_patch(mype,iopt,block_point,child_n)


! $RCSfile: amr_diagonal_patch.F90,v $
! $Revision: 1.1 $
! $Date: 2004/05/20 01:03:44 $


!------------------------------------------------------------------------
!
! This routine puts a recognizably incorrect data value into the
! extreme corner values of the solution array unk, at the beginning
! of the guard cell filling operation. At the end of the guard cell
! filling these corners are tested. The only time these values will not
! have been correctly updated is if they are diagonally opposite a
! coarser block. If this is the case then 
! a prolongation operation is needed to fill these values.
!
!
! Written :     Peter MacNeice          September 1997
!------------------------------------------------------------------------

use physicaldata
      use workspace
      use tree
      use paramesh_interfaces
      implicit none
      include  'mpif.h'


      integer mype,iopt

      real uincorrect
      parameter ( uincorrect= -1.e-30)

      real recv(nvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd)

      real recvfx(nbndvar,1:iu_bnd+1,1:ju_bnd,1:ku_bnd)
      real recvfy(nbndvar,1:iu_bnd,1:ju_bnd+k2d,1:ku_bnd)
      real recvfz(nbndvar,1:iu_bnd,1:ju_bnd,1:ku_bnd+k3d)
      real recvff(nbndvar,1:iu_bnd+1,1:ju_bnd+1,1:ku_bnd+k3d)

      integer parent_pe,parent_block
      integer ichild(2,mchild)
      integer block_point(maxblocks),child_n(2,mchild,maxblocks)
      logical lprolong

      integer imid,jmid,kmid
      integer lb
      integer jchild,ich
      integer ioff,joff,koff
      integer ia,ib,ja,jb,ka,kb
      integer nlayers

      save parent_pe,parent_block,ichild
      save recv,recvfx,recvfy,recvfz


! Apply patch-up to unk
      if(iopt.eq.1) then

! First test edge centers. Remember, in 2D edge centers become block corners.

      imid = (il_bnd+iu_bnd)/2
      jmid = (jl_bnd+ju_bnd)/2
      kmid = (kl_bnd+ku_bnd)/2


      if(lnblocks.gt.0) then
      do lb = 1,lnblocks
      if(nodetype(lb).eq.1) then

      lprolong=.false.

      if(  (unk(1,il_bnd,jl_bnd,kmid,lb).eq.uincorrect) & 
     &       .or.  (unk(1,il_bnd,ju_bnd,kmid,lb).eq.uincorrect) & 
     &       .or. (unk(1,iu_bnd,jl_bnd,kmid,lb).eq.uincorrect) & 
     &       .or.  (unk(1,iu_bnd,ju_bnd,kmid,lb).eq.uincorrect) & 
     &       ) lprolong=.true.
      if(ndim.eq.3.and.( & 
     &       (unk(1,il_bnd,jmid,kl_bnd,lb).eq.uincorrect) & 
     &       .or. (unk(1,iu_bnd,jmid,kl_bnd,lb).eq.uincorrect) & 
     &       .or. (unk(1,il_bnd,jmid,ku_bnd,lb).eq.uincorrect) & 
     &       .or. (unk(1,iu_bnd,jmid,ku_bnd,lb).eq.uincorrect) & 
     &       .or.  (unk(1,imid,jl_bnd,kl_bnd,lb).eq.uincorrect) & 
     &       .or.  (unk(1,imid,jl_bnd,ku_bnd,lb).eq.uincorrect) & 
     &       .or.  (unk(1,imid,ju_bnd,kl_bnd,lb).eq.uincorrect) & 
     &       .or.  (unk(1,imid,ju_bnd,ku_bnd,lb).eq.uincorrect) ) & 
     &       ) lprolong=.true.

      if(lprolong) then



! locate this blocks parent, copy its data and its childrens addresses
       parent_pe = parent(2,lb)
       parent_block = parent(1,lb)

       if (parent_pe.ne.mype) then
          recv(:,:,:,:) = unk(:,:,:,:,block_point(lb))
          ichild(:,:) = child_n(:,:,block_point(lb)-lnblocks)
       else
          recv(:,:,:,:) = unk(:,:,:,:,parent_block)
          ichild(:,:) = child(:,:,parent_block)
       end if

!       if(nfacevar.gt.0) then
!                call shmem_real_get
!     . (
!     .       recvfx(1,1,jl_bnd,kl_bnd),
!     .          facevarx(1,1,jl_bnd,kl_bnd,parent_block),
!     .                          len_blockfx*nbndvar,parent_pe)
!                call shmem_real_get
!     . (
!     .       recvfy(1,il_bnd,1,kl_bnd),
!     .          facevary(1,il_bnd,1,kl_bnd,parent_block),
!     .                          len_blockfy*nbndvar,parent_pe)
!                if(ndim.eq.3) then
!                call shmem_real_get
!     . (
!     .       recvfz(1,il_bnd,jl_bnd,1),
!     .          facevarz(1,il_bnd,jl_bnd,1,parent_block),
!     .                          len_blockfz*nbndvar,parent_pe)
!       endif

! identify which child this leaf block represents
       jchild = 0
       do ich=1,nchild
       if(ichild(1,ich).eq.lb.and. & 
     &       ichild(2,ich).eq.mype) jchild = ich
       enddo

! compute the offset in the parent block appropriate for this child
                ioff = mod(jchild-1,2)*nxb/2
                joff = mod((jchild-1)/2,2)*nyb/2
                koff = mod((jchild-1)/4,2)*nzb/2



      ka = 1
      kb = nzb+2*nguard*k3d

      if( unk(1,il_bnd,jl_bnd,kmid,lb).eq.uincorrect) then
       ia = 1
       ib = nguard
       ja = 1
       jb = nguard
       call amr_prolong_gen_unk_fun(recv,ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)

       if(nfacevar.gt.0) then
       recvff(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd) & 
     &       =recvfx(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
                call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb,ka,kb, & 
     &       lb,ioff,joff,koff,mype,1)

#if N_DIM >= 2
       recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd) & 
     &       =recvfy(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd)
                call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb,ka,kb, & 
     &       lb,ioff,joff,koff,mype,2)
#endif
#if N_DIM == 3
       recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d) & 
     &       =recvfz(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d)
       call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb,ka,kb+k3d, & 
     &       lb,ioff,joff,koff,mype,3)
#endif
       endif

      endif
      if( unk(1,iu_bnd,jl_bnd,kmid,lb).eq.uincorrect) then
       ia = 1+nxb+nguard
       ib = nxb+2*nguard
       ja = 1
       jb = nguard
       call amr_prolong_gen_unk_fun(recv,ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)


       if(nfacevar.gt.0) then
       recvff(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd) & 
     &       =recvfx(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
       call amr_prolong_gen_face_fun(recvff,ia+1,ib+1,ja,jb,ka,kb, & 
     &       lb,ioff,joff,koff,mype,1)

#if N_DIM >= 2
       recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd) & 
     &       =recvfy(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd)
                call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb,ka,kb, & 
     &       lb,ioff,joff,koff,mype,2)

#endif
#if N_DIM == 3
       recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d) & 
     &       =recvfz(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d)
       call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb,ka,kb+k3d, & 
     &       lb,ioff,joff,koff,mype,3)
#endif
       endif


      endif
      if( unk(1,il_bnd,ju_bnd,kmid,lb).eq.uincorrect) then
       ia =  1
       ib = nguard
       ja = 1+nyb+nguard
       jb = nyb+2*nguard
       call amr_prolong_gen_unk_fun(recv,ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)

       if(nfacevar.gt.0) then
       recvff(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd) & 
     &       =recvfx(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
                call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb,ka,kb, & 
     &       lb,ioff,joff,koff,mype,1)

#if N_DIM >= 2
       recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd) & 
     &       =recvfy(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd)
       call amr_prolong_gen_face_fun(recvff,ia,ib,ja+1,jb+1,ka,kb, & 
     &       lb,ioff,joff,koff,mype,2)
#endif
#if N_DIM == 3
       recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d) & 
     &       =recvfz(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d)
       call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb,ka,kb+k3d, & 
     &       lb,ioff,joff,koff,mype,3)
#endif
       endif



      endif
      if( unk(1,iu_bnd,ju_bnd,kmid,lb).eq.uincorrect) then
       ia = 1+nxb+nguard
       ib = nxb+2*nguard
       ja = 1+nyb+nguard
       jb = nyb+2*nguard
       call amr_prolong_gen_unk_fun(recv,ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)

       if(nfacevar.gt.0) then
       recvff(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd) & 
     &       =recvfx(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
       call amr_prolong_gen_face_fun(recvff,ia+1,ib+1,ja,jb,ka,kb, & 
     &       lb,ioff,joff,koff,mype,1)

#if N_DIM >= 2
       recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd) & 
     &       =recvfy(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd)
       call amr_prolong_gen_face_fun(recvff,ia,ib,ja+1,jb+1,ka,kb, & 
     &       lb,ioff,joff,koff,mype,2)

#endif
#if N_DIM == 3
       recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d) & 
     &       =recvfz(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d)
       call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb,ka,kb+k3d, & 
     &       lb,ioff,joff,koff,mype,3)
#endif
       endif



      endif


      if(ndim.eq.3) then

      ja = 1
      jb = nyb+2*nguard

      if( unk(1,il_bnd,jmid,kl_bnd,lb).eq.uincorrect) then
       ia = 1
       ib = nguard
       ka = 1
       kb = nguard
       call amr_prolong_gen_unk_fun(recv,ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)

       if(nfacevar.gt.0) then
       recvff(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd) & 
     &       =recvfx(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
       call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb,ka,kb, & 
     &       lb,ioff,joff,koff,mype,1)

#if N_DIM >= 2
       recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd) & 
     &       =recvfy(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd)
       call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb+1,ka,kb, & 
     &       lb,ioff,joff,koff,mype,2)

#endif
#if N_DIM == 3
       recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d) & 
     &    =recvfz(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d)
       call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb,ka,kb, & 
     &       lb,ioff,joff,koff,mype,3)
#endif
       endif



      endif
      if( unk(1,iu_bnd,jmid,kl_bnd,lb).eq.uincorrect) then
       ia = 1+nxb+nguard
       ib = nxb+2*nguard
       ka = 1
       kb = nguard
       call amr_prolong_gen_unk_fun(recv,ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)

       if(nfacevar.gt.0) then
       recvff(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd) & 
     &       =recvfx(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
       call amr_prolong_gen_face_fun(recvff,ia+1,ib+1,ja,jb,ka,kb, & 
     &       lb,ioff,joff,koff,mype,1)


#if N_DIM >= 2
       recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd) & 
     &       =recvfy(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd)
       call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb+1,ka,kb, & 
     &       lb,ioff,joff,koff,mype,2)
#endif
#if N_DIM == 3
       recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d) & 
     &       =recvfz(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d)
       call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb,ka,kb, & 
     &       lb,ioff,joff,koff,mype,3)
#endif
       endif



      endif
      if( unk(1,il_bnd,jmid,ku_bnd,lb).eq.uincorrect) then
       ia =  1
       ib = nguard
       ka = 1+nzb+nguard
       kb = nzb+2*nguard
       call amr_prolong_gen_unk_fun(recv,ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)

       if(nfacevar.gt.0) then
       recvff(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd) & 
     &       =recvfx(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
       call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb,ka,kb, & 
     &       lb,ioff,joff,koff,mype,1)

#if N_DIM >= 2
       recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd) & 
     &       =recvfy(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd)
       call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb+1,ka,kb, & 
     &       lb,ioff,joff,koff,mype,2)

#endif
#if N_DIM == 3
       recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d) & 
     &    =recvfz(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d)
       call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb,ka+k3d, & 
     &             kb+k3d,lb,ioff,joff,koff,mype,3)
#endif
       endif



      endif
      if( unk(1,iu_bnd,jmid,ku_bnd,lb).eq.uincorrect) then
       ia = 1+nxb+nguard
       ib = nxb+2*nguard
       ka = 1+nzb+nguard
       kb = nzb+2*nguard
       call amr_prolong_gen_unk_fun(recv,ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)

       if(nfacevar.gt.0) then
       recvff(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd) & 
     &       =recvfx(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
       call amr_prolong_gen_face_fun(recvff,ia+1,ib+1,ja,jb,ka,kb, & 
     &       lb,ioff,joff,koff,mype,1)
#if N_DIM >= 2

       recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd) & 
     &       =recvfy(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd)
       call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb+1,ka,kb, & 
     &       lb,ioff,joff,koff,mype,2)
#endif
#if N_DIM == 3
       recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d) & 
     &     =recvfz(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d)
       call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb,ka+k3d, & 
     &                kb+k3d,lb,ioff,joff,koff,mype,3)
#endif
       endif



      endif



      ia = 1
      ib = nxb+2*nguard

      if( unk(1,imid,jl_bnd,kl_bnd,lb).eq.uincorrect) then
       ka = 1
       kb = nguard
       ja = 1
       jb = nguard
       call amr_prolong_gen_unk_fun(recv,ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)

       if(nfacevar.gt.0) then
       recvff(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd) & 
     &       =recvfx(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
       call amr_prolong_gen_face_fun(recvff,ia,ib+1,ja,jb,ka,kb, & 
     &       lb,ioff,joff,koff,mype,1)
#if N_DIM >= 2
       recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd) & 
     &       =recvfy(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd)
       call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb,ka,kb, & 
     &       lb,ioff,joff,koff,mype,2)
#endif
#if N_DIM == 3
       recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d) & 
     &    =recvfz(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d)
       call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb,ka,kb, & 
     &       lb,ioff,joff,koff,mype,3)
#endif
       endif



      endif
      if( unk(1,imid,ju_bnd,kl_bnd,lb).eq.uincorrect) then
       ja = 1+nyb+nguard
       jb = nyb+2*nguard
       ka = 1
       kb = nguard
       call amr_prolong_gen_unk_fun(recv,ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)

       if(nfacevar.gt.0) then
       recvff(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd) & 
     &       =recvfx(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
       call amr_prolong_gen_face_fun(recvff,ia,ib+1,ja,jb,ka,kb, & 
     &       lb,ioff,joff,koff,mype,1)
#if N_DIM >= 2
       recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd) & 
     &       =recvfy(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd)
       call amr_prolong_gen_face_fun(recvff,ia,ib,ja+1,jb+1,ka,kb, & 
     &       lb,ioff,joff,koff,mype,2)
#endif
#if N_DIM == 3
       recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d) & 
     &       =recvfz(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d)
       call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb,ka,kb, & 
     &       lb,ioff,joff,koff,mype,3)
#endif
       endif



      endif
      if( unk(1,imid,jl_bnd,ku_bnd,lb).eq.uincorrect) then
       ja =  1
       jb = nguard
       ka = 1+nzb+nguard
       kb = nzb+2*nguard
       call amr_prolong_gen_unk_fun(recv,ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)

       if(nfacevar.gt.0) then
       recvff(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd) & 
     &       =recvfx(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
       call amr_prolong_gen_face_fun(recvff,ia,ib+1,ja,jb,ka,kb, & 
     &       lb,ioff,joff,koff,mype,1)
#if N_DIM >= 2
       recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd) & 
     &       =recvfy(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd)
       call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb,ka,kb, & 
     &       lb,ioff,joff,koff,mype,2)
#endif
#if N_DIM == 3
       recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d) & 
     &       =recvfz(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d)
       call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb,ka+k3d, & 
     &            kb+k3d,lb,ioff,joff,koff,mype,3)
#endif
       endif



      endif
      if( unk(1,imid,ju_bnd,ku_bnd,lb).eq.uincorrect) then
       ka = 1+nzb+nguard
       kb = nzb+2*nguard
       ja = 1+nyb+nguard
       jb = nyb+2*nguard
       call amr_prolong_gen_unk_fun(recv,ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)

       if(nfacevar.gt.0) then
       recvff(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd) & 
     &       =recvfx(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
       call amr_prolong_gen_face_fun(recvff,ia,ib+1,ja,jb,ka,kb, & 
     &       lb,ioff,joff,koff,mype,1)
#if N_DIM >= 2
       recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd) & 
     &       =recvfy(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd)
       call amr_prolong_gen_face_fun(recvff,ia,ib,ja+1,jb+1,ka,kb, & 
     &       lb,ioff,joff,koff,mype,2)
#endif
#if N_DIM == 3
       recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d) & 
     &       =recvfz(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d)
       call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb,ka+k3d, & 
     &           kb+k3d,lb,ioff,joff,koff,mype,3)
#endif
       endif



      endif


      endif ! end of ndim=3 if test



      endif ! end for if(lprolong)



      endif
      enddo
      endif


! Now test corners of blocks in 3D case.
        if(ndim.eq.3) then

        if(lnblocks.gt.0) then
        do lb = 1,lnblocks
        if(nodetype(lb).eq.1) then

        lprolong=.false.

        if(             (unk(1,il_bnd,jl_bnd,kl_bnd,lb).eq.uincorrect) & 
     &          .or.    (unk(1,il_bnd,ju_bnd,kl_bnd,lb).eq.uincorrect) & 
     &          .or.    (unk(1,iu_bnd,jl_bnd,kl_bnd,lb).eq.uincorrect) & 
     &          .or.    (unk(1,iu_bnd,ju_bnd,kl_bnd,lb).eq.uincorrect) & 
     &       .or. (unk(1,il_bnd,jl_bnd,ku_bnd,lb).eq.uincorrect) & 
     &          .or.    (unk(1,il_bnd,ju_bnd,ku_bnd,lb).eq.uincorrect) & 
     &          .or.    (unk(1,iu_bnd,jl_bnd,ku_bnd,lb).eq.uincorrect) & 
     &          .or.    (unk(1,iu_bnd,ju_bnd,ku_bnd,lb).eq.uincorrect) & 
     &          ) lprolong=.true.

        if(lprolong) then


! locate this blocks parent, copy its data and its childrens addresses
                parent_pe = parent(2,lb)
                parent_block = parent(1,lb)

                if (parent_pe.ne.mype) then
                   recv(:,:,:,:) = unk(:,:,:,:,block_point(lb))
                   ichild(:,:) = child_n(:,:,block_point(lb)-lnblocks)
                else
                   recv(:,:,:,:) = unk(:,:,:,:,parent_block)
                   ichild(:,:) = child(:,:,parent_block)
                end if

                if(nfacevar.gt.0) then
!                call shmem_real_get
!     . (
!     .       recvfx(1,1,jl_bnd,kl_bnd),
!     .          facevarx(1,1,jl_bnd,kl_bnd,parent_block),
!     .                          len_blockfx*nbndvar,parent_pe)
!                call shmem_real_get
!     . (
!     .       recvfy(1,il_bnd,1,kl_bnd),
!     .          facevary(1,il_bnd,1,kl_bnd,parent_block),
!     .                          len_blockfy*nbndvar,parent_pe)
!                call shmem_real_get
!     . (
!     .       recvfz(1,il_bnd:iu_bnd,jl_bnd:ju_bnd,1),
!     .          facevarz(1,il_bnd,jl_bnd,1,parent_block),
!     .                          len_blockfz*nbndvar,parent_pe)
                endif

! identify which child this leaf block represents
                jchild = 0
                do ich=1,nchild
                        if(ichild(1,ich).eq.lb.and. & 
     &                          ichild(2,ich).eq.mype) jchild = ich
                enddo

! compute the offset in the parent block appropriate for this child
                ioff = mod(jchild-1,2)*nxb/2
                joff = mod((jchild-1)/2,2)*nyb/2
                koff = mod((jchild-1)/4,2)*nzb/2


        ka = 1
        kb = nguard

        if( unk(1,il_bnd,jl_bnd,kl_bnd,lb).eq.uincorrect) then
                ia = 1
                ib = nguard
                ja = 1
                jb = nguard
                call amr_prolong_gen_unk_fun(recv,ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)

                if(nfacevar.gt.0) then
                recvff(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd) & 
     &          =recvfx(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
                call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb,ka,kb, & 
     &                                  lb,ioff,joff,koff,mype,1)
#if N_DIM >= 2
                recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd) & 
     &          =recvfy(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd)
                call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb,ka,kb, & 
     &                                  lb,ioff,joff,koff,mype,2)
#endif
#if N_DIM == 3
              recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d) & 
     &        =recvfz(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d)
                call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb,ka,kb, & 
     &                                  lb,ioff,joff,koff,mype,3)
#endif
                endif

        endif

        if( unk(1,il_bnd,ju_bnd,kl_bnd,lb).eq.uincorrect) then
                ia = 1
                ib = nguard
                ja = nyb+nguard+1
                jb = nyb+2*nguard
                call amr_prolong_gen_unk_fun(recv,ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)

                if(nfacevar.gt.0) then
                recvff(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd) & 
     &          =recvfx(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
                call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb,ka,kb, & 
     &                                  lb,ioff,joff,koff,mype,1)
#if N_DIM >= 2
                recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd) & 
     &          =recvfy(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd)
                call amr_prolong_gen_face_fun(recvff,ia,ib,ja+1,jb+1, & 
     &       ka,kb,lb,ioff,joff,koff,mype,2)
#endif
#if N_DIM == 3
              recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d) & 
     &        =recvfz(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d)
                call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb,ka,kb, & 
     &                                  lb,ioff,joff,koff,mype,3)
#endif
                endif

        endif

        if( unk(1,iu_bnd,jl_bnd,kl_bnd,lb).eq.uincorrect) then
                ia = nxb+nguard+1
                ib = nxb+2*nguard
                ja = 1
                jb = nguard
                call amr_prolong_gen_unk_fun(recv,ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)

                if(nfacevar.gt.0) then
                recvff(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd) & 
     &          =recvfx(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
                call amr_prolong_gen_face_fun(recvff,ia+1,ib+1,ja,jb, & 
     &       ka,kb,lb,ioff,joff,koff,mype,1)
#if N_DIM >= 2
                recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd) & 
     &          =recvfy(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd)
                call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb, & 
     &                          ka,kb,lb,ioff,joff,koff,mype,2)
#endif
#if N_DIM == 3
              recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d) & 
     &        =recvfz(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d)
                call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb,ka,kb, & 
     &                                  lb,ioff,joff,koff,mype,3)
#endif
                endif

        endif


        if( unk(1,iu_bnd,ju_bnd,kl_bnd,lb).eq.uincorrect) then
                ia = nxb+nguard+1
                ib = nxb+2*nguard
                ja = nyb+nguard+1
                jb = nyb+2*nguard
                call amr_prolong_gen_unk_fun(recv,ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)

                if(nfacevar.gt.0) then
                recvff(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd) & 
     &          =recvfx(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
                call amr_prolong_gen_face_fun(recvff,ia+1,ib+1,ja,jb, & 
     &                          ka,kb,lb,ioff,joff,koff,mype,1)
#if N_DIM >= 2
                recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd) & 
     &          =recvfy(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd)
                call amr_prolong_gen_face_fun(recvff,ia,ib,ja+1,jb+1, & 
     &                          ka,kb,lb,ioff,joff,koff,mype,2)
#endif
#if N_DIM == 3
              recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d) & 
     &        =recvfz(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d)
                call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb,ka,kb, & 
     &                                  lb,ioff,joff,koff,mype,3)
#endif
                endif

        endif

        ka = nzb+nguard+1
        kb = nzb+2*nguard

        if( unk(1,il_bnd,jl_bnd,ku_bnd,lb).eq.uincorrect) then
                ia = 1
                ib = nguard
                ja = 1
                jb = nguard
                call amr_prolong_gen_unk_fun(recv,ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)

                if(nfacevar.gt.0) then
                recvff(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd) & 
     &          =recvfx(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
                call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb,ka,kb, & 
     &                                  lb,ioff,joff,koff,mype,1)
#if N_DIM >= 2
                recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd) & 
     &          =recvfy(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd)
                call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb,ka,kb, & 
     &                                  lb,ioff,joff,koff,mype,2)
#endif
#if N_DIM == 3
              recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d) & 
     &        =recvfz(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d)
                call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb, & 
     &       ka+k3d,kb+k3d,lb,ioff,joff,koff,mype,3)
#endif
                endif

        endif

        if( unk(1,il_bnd,ju_bnd,ku_bnd,lb).eq.uincorrect) then
                ia = 1
                ib = nguard
                ja = nyb+nguard+1
                jb = nyb+2*nguard
                call amr_prolong_gen_unk_fun(recv,ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)

                if(nfacevar.gt.0) then
                recvff(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd) & 
     &          =recvfx(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
                call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb,ka,kb, & 
     &                                  lb,ioff,joff,koff,mype,1)
#if N_DIM >= 2
                recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd) & 
     &          =recvfy(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd)
                call amr_prolong_gen_face_fun(recvff,ia,ib,ja+1,jb+1, & 
     &       ka,kb,lb,ioff,joff,koff,mype,2)
#endif
#if N_DIM == 3
              recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d) & 
     &        =recvfz(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d)
                call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb, & 
     &       ka+k3d,kb+k3d,lb,ioff,joff,koff,mype,3)
#endif
                endif

        endif

        if( unk(1,iu_bnd,jl_bnd,ku_bnd,lb).eq.uincorrect) then
                ia = nxb+nguard+1
                ib = nxb+2*nguard
                ja = 1
                jb = nguard
                call amr_prolong_gen_unk_fun(recv,ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)

                if(nfacevar.gt.0) then
                recvff(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd) & 
     &          =recvfx(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
                call amr_prolong_gen_face_fun(recvff,ia+1,ib+1,ja,jb, & 
     &       ka,kb,lb,ioff,joff,koff,mype,1)
#if N_DIM >= 2
                recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd) & 
     &          =recvfy(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd)
                call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb, & 
     &                          ka,kb,lb,ioff,joff,koff,mype,2)
#endif
#if N_DIM == 3
              recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d) & 
     &        =recvfz(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d)
                call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb, & 
     &       ka+k3d,kb+k3d,lb,ioff,joff,koff,mype,3)
#endif
                endif

        endif


        if( unk(1,iu_bnd,ju_bnd,ku_bnd,lb).eq.uincorrect) then
                ia = nxb+nguard+1
                ib = nxb+2*nguard
                ja = nyb+nguard+1
                jb = nyb+2*nguard
                call amr_prolong_gen_unk_fun(recv,ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)

                if(nfacevar.gt.0) then
                recvff(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd) & 
     &          =recvfx(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
                call amr_prolong_gen_face_fun(recvff,ia+1,ib+1,ja,jb, & 
     &                          ka,kb,lb,ioff,joff,koff,mype,1)
#if N_DIM >= 2
                recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd) & 
     &          =recvfy(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd)
                call amr_prolong_gen_face_fun(recvff,ia,ib,ja+1,jb+1, & 
     &                          ka,kb,lb,ioff,joff,koff,mype,2)
#endif
#if N_DIM == 3
              recvff(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d) & 
     &        =recvfz(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d)
                call amr_prolong_gen_face_fun(recvff,ia,ib,ja,jb, & 
     &       ka+k3d,kb+k3d,lb,ioff,joff,koff,mype,3)
#endif
                endif

        endif


      endif ! end of if prolong

      endif
      enddo
      endif

      endif ! end of if ndim=3




! Apply patch-up to work
      elseif(iopt.eq.2) then

      imid = (ilw+iuw)/2
      jmid = (jlw+juw)/2
      kmid = (klw+kuw)/2

      nlayers = nguard_work


      if(lnblocks.gt.0) then
      do lb = 1,lnblocks
      if(nodetype(lb).eq.1) then

      lprolong=.false.

      if(  (work(ilw,jlw,kmid,lb,1).eq.uincorrect) & 
     &       .or.  (work(ilw,juw,kmid,lb,1).eq.uincorrect) & 
     &       .or. (work(iuw,jlw,kmid,lb,1).eq.uincorrect) & 
     &       .or.  (work(iuw,juw,kmid,lb,1).eq.uincorrect) & 
     &       ) lprolong=.true.
      if(ndim.eq.3.and.( & 
     &       (work(ilw,jmid,klw,lb,1).eq.uincorrect) & 
     &       .or. (work(iuw,jmid,klw,lb,1).eq.uincorrect) & 
     &       .or. (work(ilw,jmid,kuw,lb,1).eq.uincorrect) & 
     &       .or. (work(iuw,jmid,kuw,lb,1).eq.uincorrect) & 
     &       .or.  (work(imid,jlw,klw,lb,1).eq.uincorrect) & 
     &       .or.  (work(imid,jlw,kuw,lb,1).eq.uincorrect) & 
     &       .or.  (work(imid,juw,klw,lb,1).eq.uincorrect) & 
     &       .or.  (work(imid,juw,kuw,lb,1).eq.uincorrect) ) & 
     &       ) lprolong=.true.

      if(lprolong) then


! locate this blocks parent, copy its data and its childrens addresses
       parent_pe = parent(2,lb)
       parent_block = parent(1,lb)

       if (parent_pe.ne.mype) then
          recv1(:,:,:) = work(:,:,:,block_point(lb),1)
          ichild(:,:) = child_n(:,:,block_point(lb)-lnblocks)
       else
          recv1(:,:,:) = work(:,:,:,parent_block,1)
          ichild(:,:) = child(:,:,parent_block)
       end if

! identify which child this leaf block represents
       jchild = 0
       do ich=1,nchild
       if(ichild(1,ich).eq.lb.and. & 
     &       ichild(2,ich).eq.mype) jchild = ich
       enddo

! compute the offset in the parent block appropriate for this child
                ioff = mod(jchild-1,2)*nxb/2
                joff = mod((jchild-1)/2,2)*nyb/2
                koff = mod((jchild-1)/4,2)*nzb/2




      ka = 1
      kb = nzb+2*nguard_work*k3d

      if( work(ilw,jlw,kmid,lb,1).eq.uincorrect) then
       ia = 1
       ib = nguard_work
       ja = 1
       jb = nguard_work
       call amr_prolong_gen_work_fun(ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)


      endif
      if( work(iuw,jlw,kmid,lb,1).eq.uincorrect) then
       ia = 1+nxb+nguard_work
       ib = nxb+2*nguard_work
       ja = 1
       jb = nguard_work
       call amr_prolong_gen_work_fun(ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)


      endif
      if( work(ilw,juw,kmid,lb,1).eq.uincorrect) then
       ia =  1
       ib = nguard_work
       ja = 1+nyb+nguard_work
       jb = nyb+2*nguard_work
       call amr_prolong_gen_work_fun(ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)


      endif
      if( work(iuw,juw,kmid,lb,1).eq.uincorrect) then
       ia = 1+nxb+nguard_work
       ib = nxb+2*nguard_work
       ja = 1+nyb+nguard_work
       jb = nyb+2*nguard_work
       call amr_prolong_gen_work_fun(ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)


      endif


      if(ndim.eq.3) then

      ja = 1
      jb = nyb+2*nguard_work

      if( work(ilw,jmid,klw,lb,1).eq.uincorrect) then
       ia = 1
       ib = nguard_work
       ka = 1
       kb = nguard_work
       call amr_prolong_gen_work_fun(ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)


      endif
      if( work(iuw,jmid,klw,lb,1).eq.uincorrect) then
       ia = 1+nxb+nguard_work
       ib = nxb+2*nguard_work
       ka = 1
       kb = nguard_work
       call amr_prolong_gen_work_fun(ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)


      endif
      if( work(ilw,jmid,kuw,lb,1).eq.uincorrect) then
       ia =  1
       ib = nguard_work
       ka = 1+nzb+nguard_work
       kb = nzb+2*nguard_work
       call amr_prolong_gen_work_fun(ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)


      endif
      if( work(iuw,jmid,kuw,lb,1).eq.uincorrect) then
       ia = 1+nxb+nguard_work
       ib = nxb+2*nguard_work
       ka = 1+nzb+nguard_work
       kb = nzb+2*nguard_work
       call amr_prolong_gen_work_fun(ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)


      endif



      ia = 1
      ib = nxb+2*nguard_work

      if( work(imid,jlw,klw,lb,1).eq.uincorrect) then
       ka = 1
       kb = nguard_work
       ja = 1
       jb = nguard_work
       call amr_prolong_gen_work_fun(ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)


      endif
      if( work(imid,juw,klw,lb,1).eq.uincorrect) then
       ja = 1+nyb+nguard_work
       jb = nyb+2*nguard_work
       ka = 1
       kb = nguard_work
       call amr_prolong_gen_work_fun(ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)


      endif
      if( work(imid,jlw,kuw,lb,1).eq.uincorrect) then
       ja =  1
       jb = nguard_work
       ka = 1+nzb+nguard_work
       kb = nzb+2*nguard_work
       call amr_prolong_gen_work_fun(ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)


      endif
      if( work(imid,juw,kuw,lb,1).eq.uincorrect) then
       ka = 1+nzb+nguard_work
       kb = nzb+2*nguard_work
       ja = 1+nyb+nguard_work
       jb = nyb+2*nguard_work
       call amr_prolong_gen_work_fun(ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)


      endif


      endif ! end of ndim=3 if test



      endif ! end for if(lprolong)



      endif
      enddo ! end of loop over blocks
      endif


! Now test corners of blocks in 3D case.
        if(ndim.eq.3) then

        if(lnblocks.gt.0) then
        do lb = 1,lnblocks
        if(nodetype(lb).eq.1) then

        lprolong=.false.

        if(             (work(ilw,jlw,klw,lb,1).eq.uincorrect) & 
     &          .or.    (work(ilw,juw,klw,lb,1).eq.uincorrect) & 
     &          .or.    (work(iuw,jlw,klw,lb,1).eq.uincorrect) & 
     &          .or.    (work(iuw,juw,klw,lb,1).eq.uincorrect) & 
     &          .or.    (work(ilw,jlw,kuw,lb,1).eq.uincorrect) & 
     &          .or.    (work(ilw,juw,kuw,lb,1).eq.uincorrect) & 
     &          .or.    (work(iuw,jlw,kuw,lb,1).eq.uincorrect) & 
     &          .or.    (work(iuw,juw,kuw,lb,1).eq.uincorrect) & 
     &          ) lprolong=.true.

        if(lprolong) then


! locate this blocks parent, copy its data and its childrens addresses
                parent_pe = parent(2,lb)
                parent_block = parent(1,lb)

                if (parent_pe.ne.mype) then
                   recv1(:,:,:) = work(:,:,:,block_point(lb),1)
                   ichild(:,:) = child_n(:,:,block_point(lb)-lnblocks)
                else
                   recv1(:,:,:) = work(:,:,:,parent_block,1)
                   ichild(:,:) = child(:,:,parent_block)
                end if

! identify which child this leaf block represents
                jchild = 0
                do ich=1,nchild
                        if(ichild(1,ich).eq.lb.and. & 
     &                          ichild(2,ich).eq.mype) jchild = ich
                enddo

! compute the offset in the parent block appropriate for this child
                ioff = mod(jchild-1,2)*nxb/2
                joff = mod((jchild-1)/2,2)*nyb/2
                koff = mod((jchild-1)/4,2)*nzb/2


        ka = 1
        kb = nguard_work

        if( work(ilw,jlw,klw,lb,1).eq.uincorrect) then
                ia = 1
                ib = nguard_work
                ja = 1
                jb = nguard_work
                call amr_prolong_gen_work_fun(ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)
      endif

        if( work(ilw,juw,klw,lb,1).eq.uincorrect) then
                ia = 1
                ib = nguard_work
                ja = nyb+nguard_work+1
                jb = nyb+2*nguard_work
                call amr_prolong_gen_work_fun(ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)

        endif

        if( work(iuw,jlw,klw,lb,1).eq.uincorrect) then
                ia = nxb+nguard_work+1
                ib = nxb+2*nguard_work
                ja = 1
                jb = nguard_work
                call amr_prolong_gen_work_fun(ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)
        endif


        if( work(iuw,juw,klw,lb,1).eq.uincorrect) then
                ia = nxb+nguard_work+1
                ib = nxb+2*nguard_work
                ja = nyb+nguard_work+1
                jb = nyb+2*nguard_work
                call amr_prolong_gen_work_fun(ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)
        endif

        ka = nzb+nguard_work+1
        kb = nzb+2*nguard_work

        if( work(ilw,jlw,kuw,lb,1).eq.uincorrect) then
                ia = 1
                ib = nguard_work
                ja = 1
                jb = nguard_work
                call amr_prolong_gen_work_fun(ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)
      endif

        if( work(ilw,juw,kuw,lb,1).eq.uincorrect) then
                ia = 1
                ib = nguard_work
                ja = nyb+nguard_work+1
                jb = nyb+2*nguard_work
                call amr_prolong_gen_work_fun(ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)
        endif

        if( work(iuw,jlw,kuw,lb,1).eq.uincorrect) then
                ia = nxb+nguard_work+1
                ib = nxb+2*nguard_work
                ja = 1
                jb = nguard_work
                call amr_prolong_gen_work_fun(ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)
        endif


        if( work(iuw,juw,kuw,lb,1).eq.uincorrect) then
                ia = nxb+nguard_work+1
                ib = nxb+2*nguard_work
                ja = nyb+nguard_work+1
                jb = nyb+2*nguard_work
                call amr_prolong_gen_work_fun(ia,ib,ja,jb,ka,kb, & 
     &                          lb,ioff,joff,koff,mype)
        endif


      endif ! end of if prolong

      endif
      enddo
      endif

      endif ! end of if ndim=3






      endif ! end for if(iopt...)

      return
      end
