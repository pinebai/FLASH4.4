!#define CONSERVE
      subroutine amr_prolong_unk_fun & 
     &     (recv,isg,ioff,joff,koff,jface,mype, isrc)


! $RCSfile: amr_prolong_unk_fun.F90,v $
! $Revision: 1.2 $
! $Date: 2004/09/07 00:58:09 $


!------------------------------------------------------------------------
!
! This routine takes data from the array recv, originally extracted
! from the solution array unk, and performs a prolongation
! operation on it. The data in recv is from a parent block and the
! result of the prolongation operation is written directly into block
! isg, which is one of its children. The position of the child within the
! parent block is specified by the ioff, joff and koff arguments.
! The argument jface allows the call to limit its effect to a specific
! face of the block if required. If jface is set to a value between 1 and
! 6, then guard cells for that face are set. If jface is not between 1 to
! 6 then the prolongation operation is applied to the whole block.
!
! This particular prolongation is simple linear interpolation. It can
! be used for blocks with an even or odd number of grid cells.
!
! Conservative prolongation. Special treatment for the  cells immediately
! adjacent to a boundary (ie i=nguard,nguard+1,iu_bnd-nguard,iu_bnd-nguard+1
! and likewise for j and k indeces) if using an even number of grid cells
! per block along that axis. No special treatment is required when the number
! of cells is odd.
!
! Note: before using this routine in your program, make sure that the
! routine prolong_unk_fun_init has been called.
!
!
! Written :     Peter MacNeice          January 1997
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

      integer errorcode,ierr

!------------------------------------
! local arrays
      real recv(nvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd)

      integer isg,ioff,joff,koff,jface,mype
      integer,intent(IN),OPTIONAL :: isrc
      integer icl,icu,jcl,jcu,kcl,kcu
      integer i,j,k,ivar
      integer kp,kp1,km1,jp,jp1,jm1,ip,ip1,im1
      real wz,wy,wx,wxl,wxr,wyl,wyr,wzl,wzr

!------------------------------------

      if(prol_init.ne.100) then
       write(*,*) 'PARAMESH ERROR !'
       write(*,*) 'Error : prolong_unk_fun. ', & 
     &       'You must call amr_prolong_cc_fun_init ', & 
     &       'before you can use this routine!'
       call Driver_abortFlash("Paramesh error: amr_prolong_cc_fun_init not called");
      endif


! Set the bounds on the loop controlling the interpolation.
      icl=il_bnd
      icu=iu_bnd
      jcl=jl_bnd
      jcu=ju_bnd
      kcl=kl_bnd
      kcu=ku_bnd
      
! If interpolation is confined to 1 block face only then limit the bounds
! on the interpolation loop appropriately.
      if(jface.ge.1.and.jface.le.6) then
        icl=1+nguard
        icu=nxb+nguard
        jcl=1+nguard*k2d
        jcu=nyb+nguard*k2d
        kcl=1+nguard*k3d
        kcu=nzb+nguard*k3d
      endif
      if(jface.eq.1) then
        icl=il_bnd
        icu=il_bnd + (nguard-1)
      elseif(jface.eq.2) then
        icl=iu_bnd - (nguard-1)
        icu=iu_bnd

#if N_DIM >= 2
      elseif(jface.eq.3) then
        jcl=jl_bnd
        jcu=jl_bnd + (nguard-1)
      elseif(jface.eq.4) then
        jcl=ju_bnd - (nguard-1)
        jcu=ju_bnd
#endif
#if N_DIM == 3
      elseif(jface.eq.5) then
        kcl=kl_bnd
        kcu=kl_bnd + (nguard-1)
      elseif(jface.eq.6) then
        kcl=ku_bnd - (nguard-1)
        kcu=ku_bnd
#endif
      endif
      
      wz = 0.75-(.25**2)
      wy = wz
      wx = wz

#if N_DIM < 3
      kcl = 1
      kcu = 1
#endif


#if N_DIM < 2 
      jcl = 1
      jcu = 1
#endif


! Interpolation loop.
      do k=kcl,kcu

#if N_DIM < 3
        kp = 1
        wz = 1.
#else
        kp = ((k-1)/2)+1+(nguard/2)+koff
        kp1 = kp + k3d
        km1 = kp - k3d
#ifdef CONSERVE
! FOR CONSERVTIVE INTERPOLATION
        if (mod(nzb,2).eq.0) then
           if (k.eq.nguard.or.k.eq.nguard-1) then
              kp1 = kp
              km1 = kp
           elseif (k.eq.nguard+nzb+1.or.k.eq.nguard+nzb+2) then
              kp1 = kp
              km1 = kp
           end if
        end if
#endif
        if (mod(k,2).eq.0) then
          wzl = 0.03125               
          wzr = 0.28125
        else
          wzl = 0.28125
          wzr = 0.03125
        end if
#endif


        do j=jcl,jcu

#if N_DIM < 2
          jp = 1
          wy = 1
#else

          jp = ((j-1)/2)+1+(nguard/2)+joff
          jp1 = jp + k2d
          jm1 = jp - k2d
#ifdef CONSERVE
! FOR CONSERVTIVE INTERPOLATION
          if (mod(nyb,2).eq.0) then
             if (j.eq.nguard.or.j.eq.nguard-1) then
                jp1 = jp
                jm1 = jp
             elseif (j.eq.nguard+nyb+1.or.j.eq.nguard+nyb+2) then
                jp1 = jp
                jm1 = jp
             end if
          end if
#endif
          if (mod(j,2).eq.0) then
            wyl = 0.03125               
            wyr = 0.28125
          else
            wyl = 0.28125
            wyr = 0.03125
          end if
#endif

          do i=icl,icu
            ip = ((i-1)/2)+1+(nguard/2)+ioff
            ip1 = ip + 1
            im1 = ip - 1
#ifdef CONSERVE
            if (mod(nxb,2).eq.0) then
               if (i.eq.nguard.or.i.eq.nguard-1) then
                  ip1 = ip
                  im1 = ip
               elseif (i.eq.nguard+nxb+1.or.i.eq.nguard+nxb+2) then
                  ip1 = ip
                  im1 = ip
               end if
            end if
#endif
            if (mod(i,2).eq.0) then
              wxl = 0.03125               
              wxr = 0.28125
            else
              wxl = 0.28125
              wxr = 0.03125
            end if
! compute interpolated values at location (i,j,k)
            do ivar=1,nvar

#if N_DIM == 1
              unk(ivar,i,j,k,isg) = & 
     &           (wz*wy*wx*recv(ivar,ip,jp,kp)) + & 
     &           (wz*wy*wxl*recv(ivar,im1,jp,kp)) + & 
     &           (wz*wy*wxr*recv(ivar,ip1,jp,kp)) 
#endif

#if N_DIM == 2 
              unk(ivar,i,j,k,isg) = & 
     &           (wz*wy*wx*recv(ivar,ip,jp,kp)) + & 
     &           (wz*wy*wxl*recv(ivar,im1,jp,kp)) + & 
     &           (wz*wy*wxr*recv(ivar,ip1,jp,kp)) + &
     &           (wz*wyl*wx*recv(ivar,ip,jm1,kp)) + & 
     &           (wz*wyl*wxl*recv(ivar,im1,jm1,kp)) + & 
     &           (wz*wyl*wxr*recv(ivar,ip1,jm1,kp)) + & 
     &           (wz*wyr*wx*recv(ivar,ip,jp1,kp)) + & 
     &           (wz*wyr*wxl*recv(ivar,im1,jp1,kp)) + & 
     &           (wz*wyr*wxr*recv(ivar,ip1,jp1,kp))

#endif
#if N_DIM == 3 
              unk(ivar,i,j,k,isg) = & 
     &           (wz*wy*wx*recv(ivar,ip,jp,kp)) + & 
     &           (wz*wy*wxl*recv(ivar,im1,jp,kp)) + & 
     &           (wz*wy*wxr*recv(ivar,ip1,jp,kp)) + &
     &           (wz*wyl*wx*recv(ivar,ip,jm1,kp)) + & 
     &           (wz*wyl*wxl*recv(ivar,im1,jm1,kp)) + & 
     &           (wz*wyl*wxr*recv(ivar,ip1,jm1,kp)) + & 
     &           (wz*wyr*wx*recv(ivar,ip,jp1,kp)) + & 
     &           (wz*wyr*wxl*recv(ivar,im1,jp1,kp)) + & 
     &           (wz*wyr*wxr*recv(ivar,ip1,jp1,kp)) + &
     &           (wzl*wy*wx*recv(ivar,ip,jp,km1)) + & 
     &           (wzl*wy*wxl*recv(ivar,im1,jp,km1)) + & 
     &           (wzl*wy*wxr*recv(ivar,ip1,jp,km1)) + & 
     &           (wzl*wyl*wx*recv(ivar,ip,jm1,km1)) + & 
     &           (wzl*wyl*wxl*recv(ivar,im1,jm1,km1)) + & 
     &           (wzl*wyl*wxr*recv(ivar,ip1,jm1,km1)) + & 
     &           (wzl*wyr*wx*recv(ivar,ip,jp1,km1)) + & 
     &           (wzl*wyr*wxl*recv(ivar,im1,jp1,km1)) + & 
     &           (wzl*wyr*wxr*recv(ivar,ip1,jp1,km1)) + & 
     &           (wzr*wy*wx*recv(ivar,ip,jp,kp1)) + & 
     &           (wzr*wy*wxl*recv(ivar,im1,jp,kp1)) + & 
     &           (wzr*wy*wxr*recv(ivar,ip1,jp,kp1)) + & 
     &           (wzr*wyl*wx*recv(ivar,ip,jm1,kp1)) + & 
     &           (wzr*wyl*wxl*recv(ivar,im1,jm1,kp1)) + & 
     &           (wzr*wyl*wxr*recv(ivar,ip1,jm1,kp1)) + & 
     &           (wzr*wyr*wx*recv(ivar,ip,jp1,kp1)) + & 
     &           (wzr*wyr*wxl*recv(ivar,im1,jp1,kp1)) + & 
     &           (wzr*wyr*wxr*recv(ivar,ip1,jp1,kp1))
#endif
            enddo
          enddo
        enddo
      enddo
      
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
