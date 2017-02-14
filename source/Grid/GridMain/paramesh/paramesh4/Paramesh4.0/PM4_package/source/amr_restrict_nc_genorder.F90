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

      subroutine amr_restrict_nc_genorder(datain,dataout,ivar)




!------------------------------------------------------------------------
!
! This routine performs restriction on the array datain and
! returns the result in dataout via direct injection.
!
!
! Written :     Kevin Olson          March 2004
!------------------------------------------------------------------------


      use paramesh_dimensions
      use physicaldata

      implicit none

      real, intent(in)    :: datain(:,:,:,:)
      real, intent(inout) :: dataout(:,:,:,:)
      integer, intent(in) :: ivar

      integer :: i,j,k

!------------------------------------------------------------------------

       do k=1+nguard*k3d,nzb+nguard*k3d+k3d
       do j=1+nguard*k2d,nyb+nguard*k2d+k2d
       do i=1+nguard,nxb+nguard+1

          dataout(ivar,i,j,k) = datain(ivar,i,j,k)

       enddo
       enddo
       enddo

       return
       end subroutine amr_restrict_nc_genorder
