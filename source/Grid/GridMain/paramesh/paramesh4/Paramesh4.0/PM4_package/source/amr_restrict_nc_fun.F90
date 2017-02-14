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

      subroutine amr_restrict_nc_fun(datain,dataout)




!------------------------------------------------------------------------
!
! This routine performs restriction on the array datain and
! returns the result in dataout. Note that this does not update
! guard cell elements of dataout.
!
!
! Written :     Peter MacNeice          January 1997
!------------------------------------------------------------------------


      use paramesh_dimensions
      use physicaldata

      use paramesh_interfaces, only : amr_restrict_nc_genorder, & 
     &                                amr_restrict_nc_user

      implicit none

      real, intent(in)    :: datain(:,:,:,:)
      real, intent(inout) :: dataout(:,:,:,:)

      integer :: ivar

!------------------------------------

      if (nvarcorn > 0) then

      do ivar = 1, nvarcorn

       if (interp_mask_nc_res(ivar) < 20) then

          call amr_restrict_nc_genorder(datain,dataout,ivar)

       else

          call amr_restrict_nc_user()

       end if

      end do

      end if

      return
      end subroutine amr_restrict_nc_fun
