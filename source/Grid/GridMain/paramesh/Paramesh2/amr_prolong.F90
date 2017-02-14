      subroutine amr_prolong(mype,iopt,nlayers)


! $RCSfile: amr_prolong.F90,v $
! $Revision: 1.1 $
! $Date: 2004/05/20 01:03:44 $


!------------------------------------------------------------------------
!
! This routine interpolates data from a parent block to any
! newly created child blocks.
!
! This routine calls  routines prolong_cc and prolong_fc which
! actually perform the prolongation on the cell center and cell face
! center datastructures respectively.
!
! Arguments :
!      mype integer Current processor number
!      iopt integer Switch to select which datastructures
!       are updated. If iopt=1 the prolong_cc
!       acts on UNK, and if NBNDVAR is > 0 
!       prolong_fc updates the FACEVARX(Y)(Z)
!       arrays. If iopt=2 only WORK is updated.
!      nlayers integer number of layers of guard cells at a
!       block boundary.
!
!
! Written :     Peter MacNeice          July 1997
!------------------------------------------------------------------------
use physicaldata
      use tree
      use workspace
      implicit none
      include  'mpif.h'




#ifdef TIMINGS
#include "timer.fh"
#endif

      integer mype,iopt,nlayers

!------------------------------------

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file='amr_log',status='unknown', & 
     &        position='append')
         write (30,*) ' starting PROLONG '
         print *,' starting PROLONG '
         close (30)
      end if
#endif

      call amr_prolong_cc(mype,iopt,nlayers)

      if(nfacevar.gt.0.and.iopt.eq.1) call amr_prolong_fc(mype)

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file='amr_log',status='unknown', & 
     &        position='append')
         write (30,*) ' starting PROLONG '
         print *,' done PROLONG '
         print *,' '
         close (30)
      end if
#endif

      return
      end



