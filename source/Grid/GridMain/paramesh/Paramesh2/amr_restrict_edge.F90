      subroutine amr_restrict_edge(icoord)


! $RCSfile: amr_restrict_edge.F90,v $
! $Revision: 1.1 $
! $Date: 2004/05/20 01:03:44 $



!------------------------------------------------------------------------
!
! This routine performs a user defined reduction operation on the 
! arrays recvarx(y)(z)[1,2] and returns the result in the same arrays.
! These data arrays are defined on block boundaries only.
!
! Note that this does not update guard cell elements of recvarx(y)(z)[1,2].
!
! Also note that we use stride 2 along each dimension when computing
! reduced data values on block faces, so not all values of dataout
! have been updated.
!
!
! This particular version is only appropriate for 2nd order schemes 
! using linear interpolation with even number of mesh points along 
! each block axis.
!
! Written :     Peter MacNeice          July 1997
!------------------------------------------------------------------------
 
use physicaldata
      implicit none



!------------------------------------
      integer icoord
      integer i,j,k

!------------------------------------

      if(icoord.eq.1) then ! edges on x-face
       do k=1+nguard*k3d,nzb+(nguard+1)*k3d,2
       do j=1+nguard,nyb+nguard,2
       do i=1,2
       recvarx1e(:,i,j,k) = ( &  ! y pointing edge first
     &       recvarx1e(:,i,j,k) + & 
     &       recvarx1e(:,i,j+1,k) )*.5
       enddo
       enddo
       enddo

       if(ndim.eq.3) then
       do k=1+nguard*k3d,nzb+nguard*k3d,2
       do j=1+nguard,nyb+nguard+1,2
       do i=1,2
       recvarx2e(:,i,j,k) = ( &  ! z pointing edge 
     &       recvarx2e(:,i,j,k) + & 
     &       recvarx2e(:,i,j,k+k3d) )*.5
       enddo
       enddo
       enddo
       endif

      elseif(icoord.eq.2) then ! edges on y-face
       if(ndim.eq.3) then
       do k=1+nguard*k3d,nzb+nguard*k3d,2
       do j=1,2
       do i=1+nguard,nxb+nguard+1,2
       recvary1e(:,i,j,k) = ( &  ! z pointing edge first
     &       recvary1e(:,i,j,k) + & 
     &       recvary1e(:,i,j,k+k3d) )*.5
       enddo
       enddo
       enddo
       endif

       do k=1+nguard*k3d,nzb+(nguard+1)*k3d,2
       do j=1,2
       do i=1+nguard,nxb+nguard,2
       recvary2e(:,i,j,k) = ( &  ! x pointing edge
     &       recvary2e(:,i,j,k) + & 
     &       recvary2e(:,i+1,j,k) )*.5
       enddo
       enddo
       enddo

      elseif(icoord.eq.3) then ! edges on z-face
       do k=1,2
       do j=1+nguard,nyb+nguard+1,2
       do i=1+nguard,nxb+nguard,2
       recvarz1e(:,i,j,k) = ( &  ! x pointing edge first
     &       recvarz1e(:,i,j,k) + & 
     &       recvarz1e(:,i+1,j,k) )*.5
       enddo
       enddo
       enddo

       do k=1,2
       do j=1+nguard,nyb+nguard,2
       do i=1+nguard,nxb+nguard+1,2
       recvarz2e(:,i,j,k) = ( &  ! y pointing edge
     &       recvarz2e(:,i,j,k) + & 
     &       recvarz2e(:,i,j+1,k) )*.5
       enddo
       enddo
       enddo

      endif

      return
      end
