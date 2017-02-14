      subroutine amr_edge_average(mype,nsub)


! $RCSfile: amr_edge_average.F90,v $
! $Revision: 1.1 $
! $Date: 2004/05/20 01:03:44 $


!------------------------------------------------------------------------
! This is a wrapper routine which makes the appropriate call to the
! routines which manage edge data consistency at the boundaries between
! grid blocks of different refinement level.
!
! 
! The data structure used to store and pass this data is defined
! in the include file 'block_boundary_data.fh' which can be included
! in 'physicaldata.fh'.
!
!
! Written :     Peter MacNeice          July 1997
!------------------------------------------------------------------------

use physicaldata
      implicit none


!------------------------------------

      integer mype,nsub

#ifdef VAR_DT
      call amr_edge_average_vdt(mype,nsub) ! called if variable dt
#else
        call amr_edge_average_udt(mype)       ! called if uniform dt
      call amr_edge_diagonal_check(mype)
#endif

      return
      end
