      subroutine amr_prolong_gen_face_fun(recv,ia,ib,ja,jb,ka,kb,lb, & 
     &       ioff,joff,koff,mype,iface)


! $RCSfile: amr_prolong_gen_face_fun.F90,v $
! $Revision: 1.2 $
! $Date: 2004/09/07 00:58:09 $


!------------------------------------------------------------------------
!
! This routine takes data from the array recv, originally extracted 
! from one of the arrays facevarx(y)(z), and performs a prolongation
! operation on it. The data in recv is from a parent block and the
! result of the prolongation operation is written directly into facevarx(y)(z).
! The position of the child within the 
! parent block is specified by the ioff, joff and koff arguments.
!
! iface controls which array is updated, ie facevarx if iface=1,
! facevary if iface=2, and facevarz if iface=3.
!
! This particular prolongation is simple linear interpolation. It can
! only be used for blocks with an even number of grid cells.
!
! Conservative prolongation. Special treatment for the  cells immediately
! adjacent to a boundary (ie i=nguard,nguard+1,iu_bnd-nguard,iu_bnd-nguard+1
! and likewise for j and k indeces) if using an even number of grid cells
! per block along that axis. No special treatment is required when the number
! of cells is odd.
!
! Note: before using this routine in your program, make sure that the
! routine prolong_face_fun_init has been called.
!
!
! Written :     Peter MacNeice          July 1997
!------------------------------------------------------------------------

use physicaldata
      use tree
      implicit none
      include  'mpif.h'


      common/prolongarr_u/prol_dx(il_bnd:iu_bnd), & 
     &       prol_dy(jl_bnd:ju_bnd), prol_dz(kl_bnd:ku_bnd), & 
     &       prol_indexx(2,il_bnd:iu_bnd,2), & 
     &       prol_indexy(2,jl_bnd:ju_bnd,2), & 
     &       prol_indexz(2,kl_bnd:ku_bnd,2),prol_init
      real prol_dx,prol_dy,prol_dz
      integer prol_indexx,prol_indexy,prol_indexz
      integer prol_init

      common/prolongarr_f/prol_f_dx(il_bnd:iu_bnd+1), & 
     &       prol_f_dy(jl_bnd:ju_bnd+k2d), & 
     &       prol_f_dz(kl_bnd:ku_bnd+k3d), & 
     &       prol_f_indexx(2,il_bnd:iu_bnd+1,2), & 
     &       prol_f_indexy(2,jl_bnd:ju_bnd+k2d,2), & 
     &       prol_f_indexz(2,kl_bnd:ku_bnd+k3d,2),prol_f_init
      real prol_f_dx,prol_f_dy,prol_f_dz
      integer prol_f_indexx,prol_f_indexy,prol_f_indexz
      integer prol_f_init

!------------------------------------
! local arrays
      real recv(nbndvar,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d, & 
     &       kl_bnd:ku_bnd+k3d)

      integer errorcode,ierr

      integer ia,ib,ja,jb,ka,kb,lb,ioff,joff,koff,mype,iface
      integer icl,icu,jcl,jcu,kcl,kcu
      integer i_ind,j_ind,k_ind
      integer i,j,k,ivar
      integer i1,i1p,k1,k1p,j1,j1p
      real dx,cx,dy,cy,dz,cz

!------------------------------------

      if(prol_init.ne.100) then
       write(*,*) 'PARAMESH ERROR !'
       write(*,*) 'Error : prolong_face_fun. ', & 
     &       'You must call amr_prolong_face_fun_init ', & 
     &       'before you can use this routine!'
       call Driver_abortFlash("Paramesh error: amr_prolong_face_fun_init not called")
      endif

! Set the bounds on the loop controlling the interpolation.
        icl=ia
        icu=ib
        jcl=ja
        jcu=jb
        kcl=ka
        kcu=kb


      i_ind = 1
      j_ind = 1
      k_ind = 1
      if(ioff.gt.0) i_ind = 2
      if(joff.gt.0) j_ind = 2
      if(koff.gt.0) k_ind = 2



! Interpolation loop.
!
! Note that the range of indeces used in the facevar plane differs
! depending on the value of iface_off. This assumes that the face values
! corresponding to index 0 (ie nguard faces to the left of the block
! boundary) are never needed, when iface_off=-1. 

      if(iface.eq.1) then

        do k=kcl,kcu
             k1 = prol_indexz(1,k,k_ind)
             k1p= prol_indexz(2,k,k_ind)
             dz = prol_dz(k)
             cz = 1.-dz
             do j=jcl,jcu
                   j1 = prol_indexy(1,j,j_ind)
                   j1p= prol_indexy(2,j,j_ind)
                   dy = prol_dy(j)
                   cy = 1.-dy
                   do i=icl,icu+iface_off
                         i1 = prol_f_indexx(1,i,i_ind)
                         i1p= prol_f_indexx(2,i,i_ind)
                         dx = prol_f_dx(i)
                         cx = 1.-dx

! compute interpolated values at location (i,j,k)
                         do ivar=1,nbndvar
                             facevarx(ivar,i,j,k,lb) = & 
     &                          dz*( dy*( dx*recv(ivar,i1,j1,k1) + & 
     &                          cx*recv(ivar,i1p,j1,k1))  + & 
     &                          cy*( dx*recv(ivar,i1,j1p,k1) + & 
     &                          cx*recv(ivar,i1p,j1p,k1) ) ) + & 
     &                          cz*( dy*( dx*recv(ivar,i1,j1,k1p) + & 
     &                          cx*recv(ivar,i1p,j1,k1p))  + & 
     &                          cy*( dx*recv(ivar,i1,j1p,k1p) + & 
     &                          cx*recv(ivar,i1p,j1p,k1p) ) )
                         enddo



                    enddo
             enddo
        enddo


      elseif(iface.eq.2) then

        do k=kcl,kcu
             k1 = prol_indexz(1,k,k_ind)
             k1p= prol_indexz(2,k,k_ind)
             dz = prol_dz(k)
             cz = 1.-dz
             do j=jcl,jcu+iface_off
                   j1 = prol_f_indexy(1,j,j_ind)
                   j1p= prol_f_indexy(2,j,j_ind)
                   dy = prol_f_dy(j)
                   cy = 1.-dy
                   do i=icl,icu
                         i1 = prol_indexx(1,i,i_ind)
                         i1p= prol_indexx(2,i,i_ind)
                         dx = prol_dx(i)
                         cx = 1.-dx

! compute interpolated values at location (i,j,k)
                         do ivar=1,nbndvar
                             facevary(ivar,i,j,k,lb) = & 
     &                          dz*( dy*( dx*recv(ivar,i1,j1,k1) + & 
     &                          cx*recv(ivar,i1p,j1,k1))  + & 
     &                          cy*( dx*recv(ivar,i1,j1p,k1) + & 
     &                          cx*recv(ivar,i1p,j1p,k1) ) ) + & 
     &                          cz*( dy*( dx*recv(ivar,i1,j1,k1p) + & 
     &                          cx*recv(ivar,i1p,j1,k1p))  + & 
     &                          cy*( dx*recv(ivar,i1,j1p,k1p) + & 
     &                          cx*recv(ivar,i1p,j1p,k1p) ) )
                         enddo



                    enddo
             enddo
        enddo

      elseif(iface.eq.3) then

        do k=kcl,kcu+iface_off
             k1 = prol_f_indexz(1,k,k_ind)
             k1p= prol_f_indexz(2,k,k_ind)
             dz = prol_f_dz(k)
             cz = 1.-dz
             do j=jcl,jcu
                   j1 = prol_indexy(1,j,j_ind)
                   j1p= prol_indexy(2,j,j_ind)
                   dy = prol_dy(j)
                   cy = 1.-dy
                   do i=icl,icu
                         i1 = prol_indexx(1,i,i_ind)
                         i1p= prol_indexx(2,i,i_ind)
                         dx = prol_dx(i)
                         cx = 1.-dx

! compute interpolated values at location (i,j,k)
                         do ivar=1,nbndvar
                             facevarz(ivar,i,j,k,lb) = & 
     &                          dz*( dy*( dx*recv(ivar,i1,j1,k1) + & 
     &                          cx*recv(ivar,i1p,j1,k1))  + & 
     &                          cy*( dx*recv(ivar,i1,j1p,k1) + & 
     &                          cx*recv(ivar,i1p,j1p,k1) ) ) + & 
     &                          cz*( dy*( dx*recv(ivar,i1,j1,k1p) + & 
     &                          cx*recv(ivar,i1p,j1,k1p))  + & 
     &                          cy*( dx*recv(ivar,i1,j1p,k1p) + & 
     &                          cx*recv(ivar,i1p,j1p,k1p) ) )
                         enddo



                    enddo
             enddo
        enddo

      endif


      return
      end
