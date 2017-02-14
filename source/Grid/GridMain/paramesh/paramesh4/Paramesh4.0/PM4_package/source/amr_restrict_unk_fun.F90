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

      subroutine amr_restrict_unk_fun(datain,dataout,lb)




!------------------------------------------------------------------------
!
! This routine performs restriction on the array datain and
! returns the result in dataout. Note that this does not update
! guard cell elements of dataout.
!
! Written :     Peter MacNeice          January 1997
!------------------------------------------------------------------------


      use paramesh_dimensions
      use physicaldata

      use paramesh_interfaces, only : amr_restrict_unk_genorder, & 
     &                                amr_restrict_unk_user

      implicit none

      real, intent(in)    :: datain(:,:,:,:)
      real, intent(inout) :: dataout(:,:,:,:)
      integer, intent(in) :: lb

      integer :: ivar, order

!------------------------------------


      do ivar = 1, nvar

         if (int_gcell_on_cc(ivar)) then

         if (interp_mask_unk_res(ivar) < 20) then

! call the default interpolation routine for interpolation 
            order = interp_mask_unk_res(ivar)
            if (order <=0 .or. order > 5) order = 1
            call amr_restrict_unk_genorder(datain,dataout,order,ivar)

         elseif (interp_mask_unk_res(ivar) >= 20) then

! call a user defined routine for restriction
            call amr_restrict_unk_user()

         end if

         end if

      end do

      return
      end subroutine amr_restrict_unk_fun




