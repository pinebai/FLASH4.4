      subroutine amr_restrict_fc_fun(recv,temp,icoord)


! $RCSfile: amr_restrict_fc_fun.F90,v $
! $Revision: 1.1 $
! $Date: 2004/05/20 01:03:44 $


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
! Written :     Peter MacNeice          July 1997
!------------------------------------------------------------------------


use physicaldata
      implicit none



      real recv(nbndvar,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1, & 
     &       kl_bnd:ku_bnd+1)
      real temp(nbndvar,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1, & 
     &       kl_bnd:ku_bnd+1)

      integer icoord
      integer i,j,k,ivar


!------------------------------------

      if(icoord.eq.1) then ! x-face variables
       do k=1+nguard*k3d,nzb+nguard*k3d,2
       do j=1+nguard*k2d,nyb+nguard*k2d,2
       do i=1+nguard,nxb+nguard+1,2
       do ivar=1,nbndvar
       temp(ivar,i,j,k) = (  & 
     &       recv(ivar,i,j,k) + & 
     &       recv(ivar,i,j+k2d,k) + & 
     &       recv(ivar,i,j,k+k3d) + & 
     &       recv(ivar,i,j+k2d,k+k3d))*.25
       enddo
       enddo
       enddo
       enddo

      elseif(icoord.eq.2) then ! y-face variables
       do k=1+nguard*k3d,nzb+nguard*k3d,2
       do j=1+nguard*k2d,nyb+nguard+1,2
       do i=1+nguard,nxb+nguard,2
       do ivar=1,nbndvar
       temp(ivar,i,j,k) = (  & 
     &       recv(ivar,i,j,k) + & 
     &       recv(ivar,i+1,j,k) + & 
     &       recv(ivar,i,j,k+k3d) + & 
     &       recv(ivar,i+1,j,k+k3d))*.25
       enddo
       enddo
       enddo
       enddo

      elseif(icoord.eq.3) then ! z-face variables
       do k=1+nguard,nzb+nguard+1,2
       do j=1+nguard,nyb+nguard,2
       do i=1+nguard,nxb+nguard,2
       do ivar=1,nbndvar
       temp(ivar,i,j,k) = (  & 
     &       recv(ivar,i,j,k) + & 
     &       recv(ivar,i+1,j,k) + & 
     &       recv(ivar,i,j+1,k) + & 
     &       recv(ivar,i+1,j+1,k))*.25
       enddo
       enddo
       enddo
       enddo

      endif

      return
      end
