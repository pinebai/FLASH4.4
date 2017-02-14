      subroutine amr_mark_edges(mype,iopt)


! $RCSfile: amr_mark_edges.F90,v $
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
      implicit none
      include  'mpif.h'


      integer mype,iopt
      integer lb,jf
      integer neigh_min

      real uincorrect
      parameter ( uincorrect= -1.e-30)

      real stor1(1:iu_bnd,1:ju_bnd,1:ku_bnd)
      real storw(1:iuw,1:juw,1:kuw)


      if(iopt.eq.1) then


! identify any blocks located at an external boundary. Then zero out
! the corner values which have been set on these faces.
      if(lnblocks.gt.0) then
      do lb = 1,lnblocks
      if(nodetype(lb).eq.1) then

        neigh_min=0
        do jf=1,nfaces
          neigh_min = min(neigh_min,neigh(1,jf,lb))
        enddo


       if(neigh_min.le.-20) then

! save block face elements on boundary blocks in temporary storage array.
         stor1(il_bnd,:,:) = unk(1,il_bnd,:,:,lb)
         stor1(iu_bnd,:,:) = unk(1,iu_bnd,:,:,lb)

#if N_DIM >= 2
         stor1(:,jl_bnd,:) = unk(1,:,jl_bnd,:,lb)
         stor1(:,ju_bnd,:) = unk(1,:,ju_bnd,:,lb)
#endif
#if N_DIM == 3
           stor1(:,:,kl_bnd) = unk(1,:,:,kl_bnd,lb)
           stor1(:,:,ku_bnd) = unk(1,:,:,ku_bnd,lb)
#endif

       endif



! mark block face elements on all leaf blocks with an incorrect data value
       unk(1,il_bnd,:,:,lb) = uincorrect
       unk(1,iu_bnd,:,:,lb) = uincorrect

#if N_DIM >= 2
       unk(1,:,jl_bnd,:,lb) = uincorrect
       unk(1,:,ju_bnd,:,lb) = uincorrect
#endif
#if N_DIM == 3
         unk(1,:,:,kl_bnd,lb) = uincorrect
         unk(1,:,:,ku_bnd,lb) = uincorrect
#endif

       
       if(neigh_min.le.-20) then

!       write(*,*) 'marking: iopt1 cancel on pe/blk ',mype,lb
       if(neigh(1,1,lb).le.-20) then
         unk(1,il_bnd,:,:,lb) = stor1(1,:,:)
       endif
       if(neigh(1,2,lb).le.-20) then
         unk(1,iu_bnd,:,:,lb) = stor1(iu_bnd,:,:)
       endif

#if N_DIM >= 2
       if(neigh(1,3,lb).le.-20) then
         unk(1,:,jl_bnd,:,lb) = stor1(:,1,:)
       endif
       if(neigh(1,4,lb).le.-20) then
         unk(1,:,ju_bnd,:,lb) = stor1(:,ju_bnd,:)
       endif
#endif
#if N_DIM == 3
         if(neigh(1,5,lb).le.-20) then
           unk(1,:,:,kl_bnd,lb) = stor1(:,:,1)
         endif
         if(neigh(1,6,lb).le.-20) then
           unk(1,:,:,ku_bnd,lb) = stor1(:,:,ku_bnd)
         endif
#endif


       endif


      endif
      enddo
      endif


      elseif(iopt.eq.2) then



! identify any blocks located at an external boundary. Then zero out
! the corner values which have been set on these faces.
      if(lnblocks.gt.0) then
      do lb = 1,lnblocks
      if(nodetype(lb).eq.1) then

         neigh_min=0
         do jf=1,nfaces
           neigh_min = min(neigh_min,neigh(1,jf,lb))
         enddo


         if(neigh_min.le.-20) then
! save block face elements on boundary blocks in temporary storage array.
           storw(ilw,:,:) = work(ilw,:,:,lb,1)
           storw(iuw,:,:) = work(iuw,:,:,lb,1)

#if N_DIM >= 2
           storw(:,jlw,:) = work(:,jlw,:,lb,1)
           storw(:,juw,:) = work(:,juw,:,lb,1)
#endif
#if N_DIM == 3
             storw(:,:,klw) = work(:,:,klw,lb,1)
             storw(:,:,kuw) = work(:,:,kuw,lb,1)
#endif

         endif



! mark block face elements on all leaf blocks with an incorrect data value 
         work(ilw,:,:,lb,1) = uincorrect
         work(iuw,:,:,lb,1) = uincorrect

#if N_DIM >= 2
         work(:,jlw,:,lb,1) = uincorrect
         work(:,juw,:,lb,1) = uincorrect
#endif
#if N_DIM == 3
           work(:,:,klw,lb,1) = uincorrect
           work(:,:,kuw,lb,1) = uincorrect
#endif


       if(neigh_min.le.-20) then

!                write(*,*) 'marking: iopt2 cancel on pe/blk ',mype,lb
         if(neigh(1,1,lb).le.-20) then
           work(ilw,:,:,lb,1) = storw(ilw,:,:)
         endif
         if(neigh(1,2,lb).le.-20) then
           work(iuw,:,:,lb,1) = storw(iuw,:,:)
         endif

#if N_DIM >= 2
         if(neigh(1,3,lb).le.-20) then
           work(:,jlw,:,lb,1) = storw(:,jlw,:)
         endif
         if(neigh(1,4,lb).le.-20) then
           work(:,juw,:,lb,1) = storw(:,juw,:)
         endif
#endif
#if N_DIM == 3
           if(neigh(1,5,lb).le.-20) then
             work(:,:,klw,lb,1) = storw(:,:,klw)
           endif
           if(neigh(1,6,lb).le.-20) then
             work(:,:,kuw,lb,1) = storw(:,:,kuw)
           endif
#endif

       endif

      endif
      enddo
      endif

      endif

      return
      end
