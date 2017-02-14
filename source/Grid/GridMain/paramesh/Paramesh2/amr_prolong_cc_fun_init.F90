!#define CONSERVE       ! This selects conservative prolongation
                        ! which requires special treatment of points
                        ! immediately adjacent to the block boundaries.



      subroutine amr_prolong_cc_fun_init


! $RCSfile: amr_prolong_cc_fun_init.F90,v $
! $Revision: 1.1 $
! $Date: 2004/05/20 01:03:44 $


!------------------------------------------------------------------------
!
! This routine computes the values of dx,dy and dz used during the
! interpolation process. These are used inside the prolongation routines
! saving needless repetitive computation at the cost of minimal storage
! space.
!
! This particular prolongation is simple linear interpolation. It can
! be used for blocks with an even or odd number of grid cells.
! If CONSERVE is defined then the new mesh points immediately adjacent
! to the old block boundaries are treated specially, in a way which
! guarantees conservation.
!
! Written :     Peter MacNeice          June 1997
!------------------------------------------------------------------------

use physicaldata
      use tree
      use workspace
      implicit none
      include  'mpif.h'


! Common block to store indexing and offsets for prolongation on 'unk'.
      common/prolongarr_u/prol_dx(il_bnd:iu_bnd), & 
     &       prol_dy(jl_bnd:ju_bnd), prol_dz(kl_bnd:ku_bnd), & 
     &       prol_indexx(2,il_bnd:iu_bnd,2), & 
     &       prol_indexy(2,jl_bnd:ju_bnd,2), & 
     &          prol_indexz(2,kl_bnd:ku_bnd,2),prol_init
        real    prol_dx,prol_dy,prol_dz
        integer prol_indexx,prol_indexy,prol_indexz
        integer prol_init



! Common block to store indexing and offsets for prolongation on 'work'.
        common/prolongarr_w/prolw_dx(ilw:iuw), & 
     &          prolw_dy(jlw:juw), prolw_dz(klw:kuw), & 
     &          prolw_indexx(2,ilw:iuw,2), & 
     &          prolw_indexy(2,jlw:juw,2), & 
     &          prolw_indexz(2,klw:kuw,2),prolw_init
        real    prolw_dx,prolw_dy,prolw_dz
        integer prolw_indexx,prolw_indexy,prolw_indexz
        integer prolw_init


!------------------------------------

        real xc0,yc0,zc0
        real xc,yc,zc
        integer i,j,k,ii,jj,kk,loop
        integer jblock,jchild
        integer ioff,joff,koff
        integer i1,i1p,j1,j1p,k1,k1p

!------------------------------------

! Conditional constants. Vary depending on whether a block has an even or
! odd number of cells along a given axis.
      xc0 = .75
      yc0 = .75
      zc0 = .75
      if(mod(nxb,2).eq.1) xc0=.5
      if(mod(nyb,2).eq.1) yc0=.5
      if(mod(nzb,2).eq.1) zc0=.5


! Initialize the values of dx,dy,dz needed in prolong_unk_fun.
      do k=kl_bnd,ku_bnd
        kk = k-nguard+2*nzb
        zc = zc0+real(kk)*.5
        prol_dz(k) = mod(zc,1.)
      enddo
      do j=jl_bnd,ju_bnd
        jj = j-nguard+2*nyb
        yc = yc0+real(jj)*.5
        prol_dy(j) = mod(yc,1.)
      enddo
      do i=il_bnd,iu_bnd
        ii = i-nguard+2*nxb
        xc = xc0+real(ii)*.5
        prol_dx(i) = mod(xc,1.)
      enddo


! Compute the indeces used in the interpolation
! This includes the special conservative treatment near the block boundary.
! The outer loop selects which offset value is being used for indexing (ie
! which side of the parent block the child will be on.)

      do loop = 1,2

! compute the offset in the parent block appropriate for the different children
      if(loop.eq.1) jchild=1
      if(loop.eq.2) jchild=nchild
      ioff = mod(jchild-1,2)*nxb/2
      joff = mod((jchild-1)/2,2)*nyb/2
      koff = mod((jchild-1)/4,2)*nzb/2

                                                      ! note the 2*nxb and
      do i=il_bnd,iu_bnd                              ! nxb components in the
        ii = i-nguard                                  ! expression for i1 are
        i1 = (ii+nxb*2)/2-nxb+ioff+nguard              ! included so i1 will
        i1p= i1+1                                      ! be correct for -ve
        prol_indexx(1,i,loop) = i1                     ! values of i also.
        prol_indexx(2,i,loop) = i1p                    ! (true also for j1,k1)
#ifdef CONSERVE
        if(mod(nxb,2).eq.0) then
          if(ioff.eq.0) then
            if(i.eq.nguard) prol_indexx(2,i,loop) = i1
            if(i.eq.nguard+1) prol_indexx(1,i,loop) = i1p
          else
            if(i.eq.iu_bnd-nguard) prol_indexx(2,i,loop) = i1
            if(i.eq.iu_bnd-nguard+1) prol_indexx(1,i,loop) = i1p
          endif
        endif
#endif
      enddo

      prol_indexy(:,:,loop)=jl_bnd
      if(ndim.ge.2) then
      do j=jl_bnd,ju_bnd
        jj = j-nguard
        j1 = (jj+nyb*2)/2-nyb+joff+nguard
        j1p= j1+1
        prol_indexy(1,j,loop) = j1
        prol_indexy(2,j,loop) = j1p
#ifdef CONSERVE
        if(mod(nyb,2).eq.0) then
          if(joff.eq.0) then
            if(j.eq.nguard) prol_indexy(2,j,loop) = j1
            if(j.eq.nguard+1) prol_indexy(1,j,loop) = j1p
          else
            if(j.eq.ju_bnd-nguard) prol_indexy(2,j,loop) = j1
            if(j.eq.ju_bnd-nguard+1) prol_indexy(1,j,loop) = j1p
          endif
        endif
#endif
      enddo
      endif

      prol_indexz(:,:,loop)=kl_bnd
      if(ndim.eq.3) then
      do k=kl_bnd,ku_bnd
        kk = k-nguard
        k1 = (kk+nzb*2)/2-nzb+koff+nguard
        k1p= k1+1
        prol_indexz(1,k,loop) = k1
        prol_indexz(2,k,loop) = k1p
#ifdef CONSERVE
        if(mod(nzb,2).eq.0) then
          if(koff.eq.0) then
            if(k.eq.nguard) prol_indexz(2,k,loop) = k1
            if(k.eq.nguard+1) prol_indexz(1,k,loop) = k1p
          else
            if(k.eq.ku_bnd-nguard) prol_indexz(2,k,loop) = k1
            if(k.eq.ku_bnd-nguard+1) prol_indexz(1,k,loop) = k1p
          endif
        endif
#endif
      enddo
      endif

      enddo

! set flag to pass error check at the start of prolong_unk_fun
      prol_init = 100




! Initialize the values of dx,dy,dz needed in prolong_work_fun.
        do k=klw,kuw
          kk = k-nguard_work
          zc = zc0+real(kk)*.5
          prolw_dz(k) = mod(zc,1.)
            do j=jlw,juw
              jj = j-nguard_work
              yc = yc0+real(jj)*.5
              prolw_dy(j) = mod(yc,1.)
              do i=ilw,iuw
                ii = i-nguard_work
                xc = xc0+real(ii)*.5
                prolw_dx(i) = mod(xc,1.)
              enddo
            enddo
        enddo

! Compute the indeces used in the interpolation
! This includes the special conservative treatment near the block boundary.
! The outer loop selects which offset value is being used for indexing (ie
! which side of the parent block the child will be on.)

        do loop = 1,2

! compute the offset in the parent block appropriate for the different children
        if(loop.eq.1) jchild=1
        if(loop.eq.2) jchild=nchild
        ioff = mod(jchild-1,2)*nxb/2
        joff = mod((jchild-1)/2,2)*nyb/2
        koff = mod((jchild-1)/4,2)*nzb/2


        do i=ilw,iuw
          ii = i-nguard_work
          i1 = (ii+nxb*2)/2-nxb+ioff+nguard_work
          i1p= i1+1
          prolw_indexx(1,i,loop) = i1
          prolw_indexx(2,i,loop) = i1p
#ifdef CONSERVE
          if(mod(nxb,2).eq.0) then
            if(ioff.eq.0) then
              if(i.eq.nguard_work) prolw_indexx(2,i,loop) = i1
              if(i.eq.nguard_work+1) prolw_indexx(1,i,loop) = i1p
            else
              if(i.eq.iu_bnd-nguard_work) prolw_indexx(2,i,loop) = i1
              if(i.eq.iu_bnd-nguard_work+1) prolw_indexx(1,i,loop) = i1p
            endif
          endif
#endif
        enddo

        prolw_indexy(:,:,loop)=jlw
        if(ndim.ge.2) then
        do j=jlw,juw
          jj = j-nguard_work
          j1 = (jj+nyb*2)/2-nyb+joff+nguard_work
          j1p= j1+1
          prolw_indexy(1,j,loop) = j1
          prolw_indexy(2,j,loop) = j1p
#ifdef CONSERVE
          if(mod(nyb,2).eq.0) then
            if(joff.eq.0) then
              if(j.eq.nguard_work) prolw_indexy(2,j,loop) = j1
              if(j.eq.nguard_work+1) prolw_indexy(1,j,loop) = j1p
            else
              if(j.eq.ju_bnd-nguard_work) prolw_indexy(2,j,loop) = j1
              if(j.eq.ju_bnd-nguard_work+1) prolw_indexy(1,j,loop) = j1p
            endif
          endif
#endif
        enddo
        endif

        prolw_indexz(:,:,loop)=klw
        if(ndim.eq.3) then
        do k=klw,kuw
          kk = k-nguard_work
          k1 = (kk+nzb*2)/2-nzb+koff+nguard_work
          k1p= k1+1
          prolw_indexz(1,k,loop) = k1
          prolw_indexz(2,k,loop) = k1p
#ifdef CONSERVE
          if(mod(nzb,2).eq.0) then
            if(koff.eq.0) then
              if(k.eq.nguard_work) prolw_indexz(2,k,loop) = k1
              if(k.eq.nguard_work+1) prolw_indexz(1,k,loop) = k1p
            else
              if(k.eq.ku_bnd-nguard_work) prolw_indexz(2,k,loop) = k1
              if(k.eq.ku_bnd-nguard_work+1) prolw_indexz(1,k,loop) = k1p
            endif
          endif
#endif
        enddo
        endif

        enddo

! set flag to pass error check at the start of prolong_work_fun
      prolw_init = 100



      return
      end
