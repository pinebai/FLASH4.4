!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine amr_prolong_work_fun(ia,ib,ja,jb,ka,kb,nlayers, & 
     &       isg,ioff,joff,koff,jface,mype)


! $RCSfile: amr_prolong_work_fun.F90,v $
! $Revision: 1.2 $
! $Date: 2004/09/07 00:58:09 $


!------------------------------------------------------------------------
!
! This routine takes data from the array recv1, originally extracted
! from the workspace array work, and performs a prolongation
! operation on it, between the bounds ranges ia to ib, ja to jb, and
! ka to kb. The data in recv1 is from a parent block and the
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
! adjacent to a boundary (ie i=nguard_work,nguard_work+1,
! iu_bnd-nguard_work,iu_bnd-nguard_work+1 and likewise for j and k indeces)
! if using an even number of grid cells  per block along that axis.
! No special treatment is required when the number
! of cells is odd.
!
! Note: before using this routine in your program, make sure that the
! routine prolong_fun_init has been called.
!
!
! Written :     Peter MacNeice          January 1997
!------------------------------------------------------------------------

use physicaldata
      use tree
      use workspace
      implicit none
      include 'mpif.h'



! Common block to store indexing and offsets for prolongation on 'work'.
        common/prolongarr_w/prolw_dx(ilw:iuw), & 
     &          prolw_dy(jlw:juw), prolw_dz(klw:kuw), & 
     &          prolw_indexx(2,ilw:iuw,2), & 
     &          prolw_indexy(2,jlw:juw,2), & 
     &          prolw_indexz(2,klw:kuw,2),prolw_init
        real    prolw_dx,prolw_dy,prolw_dz
        real wx,wxl,wxr
        real wy,wyl,wyr
        real wz,wzl,wzr
        integer prolw_indexx,prolw_indexy,prolw_indexz
        integer prolw_init

        integer errorcode,ierr

        integer ia,ib,ja,jb,ka,kb,nlayers,isg,ioff,joff,koff,jface,mype
        integer icl,icu,jcl,jcu,kcl,kcu
        integer i_ind,j_ind,k_ind
        integer i,j,k
        integer ip,jp,kp

!------------------------------------

      if(prolw_init.ne.100) then
       write(*,*) 'PARAMESH ERROR !'
       write(*,*) 'Error : prolong_work_fun. ', & 
     &       'You must call amr_prolong_fun_init ', & 
     &       'before you can use this routine!'
       call Driver_abortFlash("Paramesh error: amr_prolong_fun_init not called")
      endif


! Set the bounds on the loop controlling the interpolation.
        icl=ia
        icu=ib
        jcl=ja
        jcu=jb
        kcl=ka
        kcu=kb

! If interpolation is confined to 1 block face only then limit the bounds
! on the interpolation loop appropriately.
        if(jface.ge.1.and.jface.le.6) then
             icl=1+nguard_work
             icu=nxb+nguard_work
             jcl=1+nguard_work*k2d
             jcu=nyb+nguard_work*k2d
             kcl=1+nguard_work*k3d
             kcu=nzb+nguard_work*k3d
        endif
        if(jface.eq.1) then
             icl=nguard_work+1-nlayers
             icu=nguard_work
        elseif(jface.eq.2) then
             icl=nxb+nguard_work+1
             icu=nxb+nguard_work+nlayers
        elseif(jface.eq.3) then
             jcl=nguard_work+1-nlayers
             jcu=nguard_work
        elseif(jface.eq.4) then
             jcl=nyb+nguard_work+1
             jcu=nyb+nguard_work+nlayers
        elseif(jface.eq.5) then
             kcl=nguard_work+1-nlayers
             kcu=nguard_work
        elseif(jface.eq.6) then
             kcl=nzb+nguard_work+1
             kcu=nzb+nguard_work+nlayers
        endif


      i_ind = 1
      j_ind = 1
      k_ind = 1
      if(ioff.gt.0) i_ind = 2
      if(joff.gt.0) j_ind = 2
      if(koff.gt.0) k_ind = 2

      wz = 0.75-(.25**2)
      wy = wz
      wx = wz
#if N_DIM < 3
      kcl = 1
      kcu = 1
      wz = 1.
#endif

! Interpolation loop.
        do k=kcl,kcu
#if N_DIM == 3
             kp = ((k-1)/2)+1+(nguard_work/2)+koff
             if (mod(k,2).eq.0) then
               wzl = 0.03125               
               wzr = 0.28125
             else
               wzl = 0.28125
               wzr = 0.03125
             end if
#else
             kp = 1
#endif
             do j=jcl,jcu
               jp = ((j-1)/2)+1+(nguard_work/2)+joff
               if (mod(j,2).eq.0) then
                 wyl = 0.03125               
                 wyr = 0.28125
               else
                 wyl = 0.28125
                 wyr = 0.03125
               end if
               do i=icl,icu
                 ip = ((i-1)/2)+1+(nguard_work/2)+ioff
                 if (mod(i,2).eq.0) then
                   wxl = 0.03125               
                   wxr = 0.28125
                 else
                   wxl = 0.28125
                   wxr = 0.03125
                 end if
! compute interpolated values at location (i,j,k)
#if N_DIM < 3
                 work(i,j,k,isg,1) = & 
     &                          (wz*wy*wx*recv1(ip,jp,kp)) + & 
     &                          (wz*wy*wxl*recv1(ip-1,jp,kp)) + & 
     &                          (wz*wy*wxr*recv1(ip+1,jp,kp)) + & 
     &                          (wz*wyl*wx*recv1(ip,jp-1,kp)) + & 
     &                          (wz*wyl*wxl*recv1(ip-1,jp-1,kp)) + & 
     &                          (wz*wyl*wxr*recv1(ip+1,jp-1,kp)) + & 
     &                          (wz*wyr*wx*recv1(ip,jp+1,kp)) + & 
     &                          (wz*wyr*wxl*recv1(ip-1,jp+1,kp)) + & 
     &                          (wz*wyr*wxr*recv1(ip+1,jp+1,kp))
#endif
#if N_DIM == 3  
                 work(i,j,k,isg,1) = & 
     &                          (wz*wy*wx*recv1(ip,jp,kp)) + & 
     &                          (wz*wy*wxl*recv1(ip-1,jp,kp)) + & 
     &                          (wz*wy*wxr*recv1(ip+1,jp,kp)) + & 
     &                          (wz*wyl*wx*recv1(ip,jp-1,kp)) + & 
     &                          (wz*wyl*wxl*recv1(ip-1,jp-1,kp)) + & 
     &                          (wz*wyl*wxr*recv1(ip+1,jp-1,kp)) + & 
     &                          (wz*wyr*wx*recv1(ip,jp+1,kp)) + & 
     &                          (wz*wyr*wxl*recv1(ip-1,jp+1,kp)) + & 
     &                          (wz*wyr*wxr*recv1(ip+1,jp+1,kp)) + &
     &                          (wzl*wy*wx*recv1(ip,jp,kp-1)) + & 
     &                          (wzl*wy*wxl*recv1(ip-1,jp,kp-1)) + & 
     &                          (wzl*wy*wxr*recv1(ip+1,jp,kp-1)) + & 
     &                          (wzl*wyl*wx*recv1(ip,jp-1,kp-1)) + & 
     &                          (wzl*wyl*wxl*recv1(ip-1,jp-1,kp-1)) + & 
     &                          (wzl*wyl*wxr*recv1(ip+1,jp-1,kp-1)) + & 
     &                          (wzl*wyr*wx*recv1(ip,jp+1,kp-1)) + & 
     &                          (wzl*wyr*wxl*recv1(ip-1,jp+1,kp-1)) + & 
     &                          (wzl*wyr*wxr*recv1(ip+1,jp+1,kp-1)) + & 
     &                          (wzr*wy*wx*recv1(ip,jp,kp+1)) + & 
     &                          (wzr*wy*wxl*recv1(ip-1,jp,kp+1)) + & 
     &                          (wzr*wy*wxr*recv1(ip+1,jp,kp+1)) + & 
     &                          (wzr*wyl*wx*recv1(ip,jp-1,kp+1)) + & 
     &                          (wzr*wyl*wxl*recv1(ip-1,jp-1,kp+1)) + & 
     &                          (wzr*wyl*wxr*recv1(ip+1,jp-1,kp+1)) + & 
     &                          (wzr*wyr*wx*recv1(ip,jp+1,kp+1)) + & 
     &                          (wzr*wyr*wxl*recv1(ip-1,jp+1,kp+1)) + & 
     &                          (wzr*wyr*wxr*recv1(ip+1,jp+1,kp+1))
#endif

                    enddo
             enddo
        enddo

      return
      end

