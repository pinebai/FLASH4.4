      subroutine amr_prolong_gen_work_fun(ia,ib,ja,jb,ka,kb, & 
     &       lb,ioff,joff,koff,mype)


! $RCSfile: amr_prolong_gen_work_fun.F90,v $
! $Revision: 1.2 $
! $Date: 2004/09/07 00:58:09 $


!------------------------------------------------------------------------
!
! This routine takes data from the array recv, originally extracted 
! from the workspace array work on some block, 
! and performs a prolongation operation on it, between the bounds ranges 
! ia to ib, ja to jb, and ka to kb. The data in recv1 is from a parent 
! block and the result of the prolongation operation is returned in
! outp.
! The position of the child within the 
! parent block is specified by the ioff, joff and koff arguments.
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
        integer prolw_indexx,prolw_indexy,prolw_indexz
        integer prolw_init

        integer errorcode,ierr
        
        integer ia,ib,ja,jb,ka,kb,lb,ioff,joff,koff,mype
        integer icl,icu,jcl,jcu,kcl,kcu
        integer i_ind,j_ind,k_ind
        integer i1,i1p,j1,j1p,k1,k1p
        integer i,j,k

        real dx,cx,dy,cy,dz,cz

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


      i_ind = 1
      j_ind = 1
      k_ind = 1
      if(ioff.gt.0) i_ind = 2
      if(joff.gt.0) j_ind = 2
      if(koff.gt.0) k_ind = 2

! Interpolation loop.
        do k=kcl,kcu
             k1 = prolw_indexz(1,k,k_ind)
             k1p= prolw_indexz(2,k,k_ind)
             dz = prolw_dz(k)
             cz = 1.-dz
             do j=jcl,jcu
                   j1 = prolw_indexy(1,j,j_ind)
                   j1p= prolw_indexy(2,j,j_ind)
                   dy = prolw_dy(j)
                   cy = 1.-dy
                   do i=icl,icu
                         i1 = prolw_indexx(1,i,i_ind)
                         i1p= prolw_indexx(2,i,i_ind)
                         dx = prolw_dx(i)
                         cx = 1.-dx

! compute interpolated values at location (i,j,k)
                             work(i,j,k,lb,1) = & 
     &                          dz*( dy*( dx*recv1(i1,j1,k1) + & 
     &                          cx*recv1(i1p,j1,k1))  + & 
     &                          cy*( dx*recv1(i1,j1p,k1) + & 
     &                          cx*recv1(i1p,j1p,k1) ) ) + & 
     &                          cz*( dy*( dx*recv1(i1,j1,k1p) + & 
     &                          cx*recv1(i1p,j1,k1p))  + & 
     &                          cy*( dx*recv1(i1,j1p,k1p) + & 
     &                          cx*recv1(i1p,j1p,k1p) ) )



                    enddo
             enddo
        enddo


      return
      end
