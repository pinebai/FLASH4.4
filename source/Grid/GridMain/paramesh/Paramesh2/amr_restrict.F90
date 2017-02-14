      subroutine amr_restrict(mype,iopt,iempty)


! $RCSfile: amr_restrict.F90,v $
! $Revision: 1.1 $
! $Date: 2004/05/20 01:03:44 $


!------------------------------------------------------------------------
!
! This routine does the data averaging required when a child block
! passes data back to its parent. The parent receives interior data
! only, not guard cell data.
! This routine calls a user provided routine called restrict_fun
! which defines the pattern of restriction which the user wishes to
! apply.
!
! Arguments :
!       mype    integer         Current processor number
!       iopt    integer         Switch to select which datastructures
!                                are updated. If iopt=1 the restrict_cc
!                                acts on UNK, and if NBNDVAR is > 0 
!                                restrict_fc updates the FACEVARX(Y)(Z)
!                                arrays. If iopt=2 only WORK is updated.
!       iempty  integer         Switch to choose between unlimited 
!                                scope of the restriction operation(=0)
!                                and limiting the restriction to deal
!                                only with leaf blocks bordering on
!                                coarser empty leaf blocks(=1)
!
!
! Written :     Peter MacNeice          July 1997
!------------------------------------------------------------------------

use physicaldata
      use tree
      implicit none
      include 'mpif.h'



#ifdef TIMINGS
#include "timer.fh"
#endif

      integer, intent(in) :: mype,iopt,iempty
      integer i

! local arrays
        logical rflag(maxblocks)
      save rflag
!------------------------------------------------------------------------

        do i = 1, lnblocks
          rflag(i) = .false.
        enddo


! If this restriction call amr_is to prepare to deal with empty grid
! cells then limit the scope of the restriction to leaf blocks
! bordering coarser empty blocks.
        if(iempty.eq.1) then
                call amr_restrict_eblock_marker(rflag)
        else
! otherwise allow unlimited restriction
                do i = 1, lnblocks
                  rflag(i) = .true.
                enddo
        endif


! Now perform selected restrictions.

      call amr_restrict_cc(mype,iopt,rflag)

      if(nfacevar.gt.0.and.iopt.eq.1) call amr_restrict_fc(mype,rflag)


      return
      end


