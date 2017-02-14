      subroutine amr_prolong_fun_init


! $RCSfile: amr_prolong_fun_init.F90,v $
! $Revision: 1.1 $
! $Date: 2004/05/20 01:03:44 $

      implicit none

      include 'mpif.h'


!------------------------------------------------------------------------
!
! This routine calls the routines which compute the values of dx,dy and 
! dz and some index vectors used during the interpolation process. 
! These are used inside the prolongation routines
! saving needless repetitive computation at the cost of minimal storage
! space.
!
! Written :     Peter MacNeice          July 1997
!------------------------------------------------------------------------

      call amr_prolong_cc_fun_init

      call amr_prolong_face_fun_init


      return
      end

