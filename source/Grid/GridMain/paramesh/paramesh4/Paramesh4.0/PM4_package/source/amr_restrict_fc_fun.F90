!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"


      subroutine amr_restrict_fc_fun(recv,temp,icoord)




!------------------------------------------------------------------------
!
! This routine performs a user defined reduction operation on the 
! array recv and returns the result in temp.
!
! Note that this does not update guard cell elements of temp.
!
! Also note that we use stride 2 along each dimension when computing
! reduced data values on cell faces, so not all values of temp
! have been updated.
!
!
! Written :     Peter MacNeice          July 1997
!!
!!    Modified to call amr_restrict_fc_mg   Klaus Weide  Sep 2016
!------------------------------------------------------------------------


      use paramesh_dimensions
      use physicaldata

      use paramesh_interfaces, only : amr_restrict_fc_genorder, & 
                                      amr_restrict_fc_ins,             &
                                      amr_restrict_fc_user

      implicit none

      real,    intent(in)    :: recv(:,:,:,:)
      real,    intent(inout) :: temp(:,:,:,:)
      integer, intent(in)    :: icoord

      integer :: order, ivar


      do ivar = 1, nbndvar

         if (icoord == 1) then
            order = interp_mask_facex_res(ivar)
         elseif (icoord == 2) then
            order = interp_mask_facey_res(ivar)
         elseif (icoord == 3) then
            order = interp_mask_facez_res(ivar)
         end if

         if (order < 20) then

            if (order <= 0 .or. order > 5) order = 1
            call amr_restrict_fc_genorder(recv,temp,icoord,order,ivar)

         else if (order .GE. 30 .AND. order .LE. 35) then

            call amr_restrict_fc_ins(recv,temp,icoord,order-30,ivar)

         else

            call amr_restrict_fc_user()

         end if

      end do

      return
      end subroutine amr_restrict_fc_fun

