      subroutine amr_prolong_gen_unk_fun(recv,ia,ib,ja,jb,ka,kb, & 
     &       isg,ioff,joff,koff,mype, isrc)


! $RCSfile: amr_prolong_gen_unk_fun.F90,v $
! $Revision: 1.2 $
! $Date: 2004/09/07 00:58:09 $


! 2ND order interpolation !!!
        
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

!------------------------------------
! local arrays
      real recv(nvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd)

      integer ia,ib,ja,jb,ka,kb,isg,ioff,joff,koff,mype
      integer,intent(IN),OPTIONAL :: isrc

      integer icl,icu,jcl,jcu,kcl,kcu
      integer i,j,k,ivar
      integer kp,kp1,jp,jp1,ip,ip1,im1,jm1,km1
      real wx,wy,wz
      real wzl,wzr,wyl,wyr,wxl,wxr

      integer errorcode,ierr


!------------------------------------

      if(prol_init.ne.100) then
       write(*,*) 'PARAMESH ERROR !'
       write(*,*) 'Error : prolong_gen_unk_fun. ', & 
     &       'You must call amr_initialize_amr ', & 
     &       'before you can use this routine!'
       call Driver_abortFlash("Paramesh error: amr_initialize_amr not called");
      endif


! Set the bounds on the loop controlling the interpolation.
        icl=ia
        icu=ib
        jcl=ja
        jcu=jb
        kcl=ka
        kcu=kb

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
     &          (wzl*wy*wx*recv(ivar,ip,jp,km1)) + & 
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
