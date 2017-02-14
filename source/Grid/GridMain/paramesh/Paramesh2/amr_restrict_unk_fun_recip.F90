      subroutine amr_restrict_unk_fun_recip(datain,dataout)


! $RCSfile: amr_restrict_unk_fun_recip.F90,v $
! $Revision: 1.1 $
! $Date: 2004/05/20 01:03:44 $


!------------------------------------------------------------------------
!
! This routine performs restriction on the array datain and
! returns the result in dataout. Note that this does not update
! guard cell elements of dataout.
!
! This particular version applies the 3D generalization of the
! restriction operator in eqn (19.6.17) of the 2nd edition of
! Numerical recipes.
! The 2D case is
!                  | 1  2  1 |
!                  | 2  4  2 | /16.
!                  | 1  2  1 |
!
!
! Written :     Peter MacNeice          January 1997
!------------------------------------------------------------------------

use physicaldata
      implicit none


!------------------------------------
! local arrays
      real datain(nvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
      real dataout(nvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd)

      integer i,j,k,ivar

!------------------------------------

       do k=1+nguard*k3d,nzb+nguard*k3d
       do j=1+nguard*k2d,nyb+nguard*k2d
       do i=1+nguard,nxb+nguard
       do ivar=1,nvar
       dataout(ivar,i,j,k) = ( & 
     &   ( datain(ivar,i-1,j-k2d,k-k3d) + 2.*datain(ivar,i,j-k2d,k-k3d) & 
     &              + datain(ivar,i+1,j-k2d,k-k3d) ) + & 
     &  2.*( datain(ivar,i-1,j,k-k3d) + 2.*datain(ivar,i,j,k-k3d) + & 
     &              datain(ivar,i+1,j,k-k3d) ) + & 
     &   ( datain(ivar,i-1,j+k2d,k-k3d) + 2.*datain(ivar,i,j+k2d,k-k3d) & 
     &              + datain(ivar,i+1,j+k2d,k-k3d) ) + & 
     &  2.*( datain(ivar,i-1,j-k2d,k) + 2.*datain(ivar,i,j-k2d,k) + & 
     &              datain(ivar,i+1,j-k2d,k) ) + & 
     &  4.*( datain(ivar,i-1,j,k) + 2.*datain(ivar,i,j,k) + & 
     &              datain(ivar,i+1,j,k) ) + & 
     &  2.*( datain(ivar,i-1,j+k2d,k) + 2.*datain(ivar,i,j+k2d,k) + & 
     &              datain(ivar,i+1,j+k2d,k) ) + & 
     &   ( datain(ivar,i-1,j-k2d,k+k3d) + 2.*datain(ivar,i,j-k2d,k+k3d) & 
     &              + datain(ivar,i+1,j-k2d,k+k3d) ) + & 
     &  2.*( datain(ivar,i-1,j,k+k3d) + 2.*datain(ivar,i,j,k+k3d) + & 
     &              datain(ivar,i+1,j,k+k3d) ) + & 
     &   ( datain(ivar,i-1,j+k2d,k+k3d) + 2.*datain(ivar,i,j+k2d,k+k3d) & 
     &           +   datain(ivar,i+1,j+k2d,k+k3d) )      )/64.

       enddo
       enddo
       enddo
       enddo

      return
      end
