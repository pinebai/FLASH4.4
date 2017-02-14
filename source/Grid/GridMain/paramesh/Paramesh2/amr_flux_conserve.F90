      subroutine amr_flux_conserve(mype,nsub,idir)


! $RCSfile: amr_flux_conserve.F90,v $
! $Revision: 1.1 $
! $Date: 2004/05/20 01:03:44 $


!------------------------------------------------------------------------
! This is a wrapper routine which makes the appropriate call to the
! routines which manage flux conservation at the boundaries between
! grid blocks of different refinement level.
!
! 
! These routines get block boundary data from neighbors who are
! parents of leaf blocks. This is required in flux conserving schemes
! where the coarser block needs to use the same fluxes and mean pressures
! as will be used on the finer blocks across their shared boundary.
!
! The data structure used to store and pass this data is defined
! in the include file 'block_boundary_data.h' which can be included
! in 'physicaldata.h'.
!
!
! Written :     Peter MacNeice          February 1997
!------------------------------------------------------------------------

use physicaldata
      use tree
      implicit none
      include 'mpif.h'


!------------------------------------

      integer mype,nsub,idir

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file='amr_log',status='unknown', & 
     &        position='append')
         write (30,*) ' starting FLUX_CONSERVE '
         print *,' starting FLUX_CONSERVE '
         close (30)
      end if
#endif

#ifdef VAR_DT
      call amr_flux_conserve_vdt(mype,nsub) ! called if variable dt
#else
      call amr_flux_conserve_udt(mype,idir) ! called if uniform dt
#endif

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file='amr_log',status='unknown', & 
     &        position='append')
         write (30,*) ' done FLUX_CONSERVE '
         write (30,*) ' '
         print *,' done FLUX_CONSERVE '
         print *,' '
         close (30)
      end if
#endif

      return
      end
