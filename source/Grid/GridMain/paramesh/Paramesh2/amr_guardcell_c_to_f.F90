      subroutine amr_guardcell_c_to_f(mype,iopt,nlayers,idiag,idir)


! $RCSfile: amr_guardcell_c_to_f.F90,v $
! $Revision: 1.1 $
! $Date: 2004/05/20 01:03:44 $


!------------------------------------------------------------------------
! If no cell-face-centered data is being used then this call
! will work for both odd or even numbers of cells along any
! block coordinate axis. However if cell-face-centered data is 
! being used then it will only work for even block sizes.
!
! This routine manages the transfer of guard cell data to all
! leaf blocks from any neighbors which they have at a coarser 
! resolution.
!
! Written :     Peter MacNeice          July 1997
!------------------------------------------------------------------------
!
! Arguments:
!       mype            local processor number
!       iopt            a switch to control which data source is to be used
!                       iopt=1 will use 'unk'
!                       iopt=2 will use 'work'
!       nlayers         the number of guard cell layers at each boundary
!
!       idiag           get diagonal guardcells if idiag = 1
!
!------------------------------------

use physicaldata
      use tree
      implicit none
      include 'mpif.h'


      integer mype,iopt,nlayers,idiag,idir

      integer block_point(maxblocks),child_n(2,mchild,maxblocks)


      call amr_guardcell_cc_c_to_f(mype,iopt,nlayers,idir, & 
     &     block_point,child_n)

      if(nfacevar.gt.0.and.iopt.eq.1)  & 
     &       call amr_guardcell_fc_c_to_f(mype)


      if (idiag.eq.1) then
        call amr_diagonal_patch(mype,iopt,block_point,child_n)
      end if

      return
      end
