!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


      subroutine amr_1blk_cc_prol_gen_work_fun(recvt, & 
     &       ia,ib,ja,jb,ka,kb, & 
     &       idest,ioff,joff,koff,mype,lb,pe_p,lb_p,interp)



!
!------------------------------------------------------------------------
!
! This routine is a wrapper routine which calls the functions
! which prolong data for WORK. The argument idest selects the layer
! within WORK on which prolongation is required.
! The argument interp selects the actual interpolation function to be
! called.
!
!
! Written :     Peter MacNeice          January 2002
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace
      use prolong_arrays

      use paramesh_interfaces, only :  & 
     &                  amr_1blk_cc_prol_work_inject, & 
     &                  amr_1blk_cc_prol_work_linear, & 
     &                  amr_1blk_cc_prol_work_user, & 
     &                  amr_1blk_cc_prol_work_genorder,&
                        amr_1blk_cc_prol_work_mg


      implicit none

      integer, intent(in) :: ia,ib,ja,jb,ka,kb,idest
      integer, intent(in) :: ioff,joff,koff,mype
      integer, intent(in) :: lb,lb_p,pe_p
      integer, intent(in) :: interp
      real,    intent(inout) :: recvt(:,:,:)

!------------------------------------

      if (interp < 20) then

! Simple Injection
      if(interp.eq.0) & 
     &  call amr_1blk_cc_prol_work_inject(recvt, & 
     &       ia,ib,ja,jb,ka,kb, & 
     &       idest,ioff,joff,koff,mype)

! Linear interpolation
      if(interp.eq.1) & 
     &   call amr_1blk_cc_prol_work_linear(recvt, & 
     &       ia,ib,ja,jb,ka,kb, & 
     &       idest,ioff,joff,koff,mype)

! General Order
      
      if (interp > 1)  & 
     &   call amr_1blk_cc_prol_work_genorder(recvt, & 
     &       ia,ib,ja,jb,ka,kb, & 
     &       idest,ioff,joff,koff,mype,interp)

      elseif (interp == 20) then

      call amr_1blk_cc_prol_work_user()

      elseif (interp .GE. 30 .AND. interp .LE. 32) then

         call amr_1blk_cc_prol_work_mg(recvt,       &
             ia,ib,ja,jb,ka,kb,                     &
             idest,ioff,joff,koff,lb,interp-30)

      end if

      return
      end subroutine amr_1blk_cc_prol_gen_work_fun
