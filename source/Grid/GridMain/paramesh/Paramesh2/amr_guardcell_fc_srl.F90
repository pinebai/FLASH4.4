      subroutine amr_guardcell_fc_srl(mype,nlayers)


! $RCSfile: amr_guardcell_fc_srl.F90,v $
! $Revision: 1.1 $
! $Date: 2004/05/20 01:03:44 $



!------------------------------------------------------------------------
!
! This routine manages the exchange of guard cell information between
! blocks assuming that exchange is only required between blocks at the
! same refinement level. This routine operates on idata defined at
! cell face centers (ie in facevarx etc).
! The actual exchanges are performed with calls to the routines 
! face_cp_loc and face_cp_remote.
!
! The order of the loops (ie placing the loop over faces on the outside,
! with a shmem_barrier_all after each iteration) is required to guarantee that the
! guardcells along block edges are correct.
!
!
! Written :     Peter MacNeice          July 1997
!------------------------------------------------------------------------
!
! Arguments:
!      mype local processor number
!      nlayers the number of guard cell layers at each boundary
!
!------------------------------------

use physicaldata
      use tree
      use workspace
      implicit none


!------------------------------------

      integer mype,nlayers

      integer j,isg
      integer jneigh,jblock

! cycle through block faces
      do j = 1,nfaces

! cycle through the grid blocks on this processor
      if(lnblocks.gt.0) then
      do isg = 1,lnblocks

! if current block has no children or if it has at least one leaf
! child then get the guard cell info that it needs.
      if(nodetype(isg).le.2) then

       jneigh = neigh(1,j,isg)

! If neighbor block exists at this refinement level
       if(jneigh.gt.-1) then

       if(neigh(2,j,isg).eq.mype) then
       call amr_face_cp_loc(isg,jneigh,j, & 
     &       nlayers)
       else
       call amr_face_cp_remote(mype,neigh(2,j,isg), & 
     &       jneigh,isg,j,nlayers)
       endif
       endif


      endif

      enddo
      endif

!       if(j.eq.2.or.j.eq.4) call shmem_barrier_all()
!      write(*,*) 'proc ',mype,' cleared j= ',j,' shmem_barrier_all'
      enddo

!       call shmem_barrier_all()
!      write(*,*) 'proc ',mype,' cleared j= 6 shmem_barrier_all'

      return
      end
