      subroutine amr_guardcell_srl(mype,iopt,nlayers,idiag,idir,maxNodetype_gcWanted)


! $RCSfile: amr_guardcell_srl.F90,v $
! $Revision: 1.1 $
! $Date: 2004/05/20 01:03:44 $



!------------------------------------------------------------------------
!
! This routine manages the exchange of guard cell information between
! blocks assuming that exchange is only required between blocks at the
! same refinement level. Calls are made to routines guardcell_cc_srl
! and guardcell_cc_srl which perform the data movement for cell centered
! and face centered data respectively.
!
! For general guard cell filling use iopt=1. This means the cell centered
! array unk and the face centered arrays facevarx(y)(z) are filled. If
! iopt=2 then only the cell centered array work is filled.
!
! Written :     Peter MacNeice          July 1997
!------------------------------------------------------------------------
!
! Arguments:
!      mype local processor number
!      iopt a switch to control which data source is to be used
!       iopt=1 will use 'unk'
!       iopt=2 will use 'work'
!      nlayers the number of guard cell layers at each boundary
!
!------------------------------------

use physicaldata
      use tree
      use workspace
      implicit none
      include  'mpif.h'




#ifdef TIMINGS
#include "timer.fh"
#endif
      integer mype,iopt,nlayers,idiag,idir
      integer,intent(IN) :: maxNodetype_gcWanted

!------------------------------------

      call amr_guardcell_cc_srl(mype,iopt,nlayers,idiag,idir,maxNodetype_gcWanted)


      if(nfacevar.gt.0.and.iopt.eq.1)  & 
     &       call amr_guardcell_fc_srl(mype,nlayers)


      return
      end
