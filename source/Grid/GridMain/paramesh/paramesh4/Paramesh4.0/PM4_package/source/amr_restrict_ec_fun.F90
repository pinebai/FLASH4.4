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


      subroutine amr_restrict_ec_fun(recv,temp,icoord)




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
! This particular version is only appropriate for 2nd order schemes 
! using linear interpolation with even number of mesh points along 
! each block axis.
!
! Written :     Peter MacNeice          December 2000
!------------------------------------------------------------------------


      use paramesh_dimensions
      use physicaldata

      use paramesh_interfaces, only : amr_restrict_ec_genorder, & 
     &                                amr_restrict_ec_user

      implicit none

      integer, intent(in)    :: icoord
      real,    intent(in)    :: recv(:,:,:,:)
      real,    intent(inout) :: temp(:,:,:,:)

      integer :: ivar, order



!------------------------------------


      do ivar = 1, nbndvare
         
         order = interp_mask_ec_res(ivar)
         
         if (order < 20) then

            if (order <= 0 .or. order > 5) order = 1
            call amr_restrict_ec_genorder(recv,temp,icoord,order,ivar)
            
         else
            
            call amr_restrict_ec_user()
            
         end if
         
      end do

      return
      end subroutine amr_restrict_ec_fun
