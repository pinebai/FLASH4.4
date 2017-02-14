!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/mpi_amr_1blk_guardcell_c_to_f
!! NAME
!!
!!   mpi_amr_1blk_guardcell_c_to_f
!!
!! SYNOPSIS
!!
!!   call mpi_amr_1blk_guarcell_c_to_f(mype, lb, pe, iopt, nlayers, surrblks,
!!                                     lcc,lfc,lec,lnc,icoord,ldiag,
!!                                     nlayersx,nlayersy,nlayersz) 
!!
!!   call mpi_amr_1blk_guarcell_c_to_f(integer, integer, integer, integer, integer, 
!!                                     integer,
!!                                     logical, logical, logical, logical, integer, 
!!                                     logical,
!!                                     integer, integer, integer) 
!!
!! ARGUMENTS
!!
!!   integer, intent(in) :: mype           
!!        The local processor number.
!!
!!   integer, intent(in) :: lb
!!        The selected block
!!
!!   integer, intent(in) :: pe
!!        Processor storing the selected block
!!
!!   integer, intent(in) :: iopt     
!!        A switch to control which data source is to be used, iopt = 1 will use 'unk',
!!        'facevarx,y,z', 'unk_e_x,y,z', and 'unk_n'.  iopt >= 2 will use 'work'.
!!
!!   integer, intent(in) :: nlayers
!!        The number of guard cell layers at each boundary.  NOTE: This no longer has
!!        any effect.  Use nlayersx, y and z (see below).
!!
!!   integer, intent(in) :: surrblks(:,:,:,:)
!!        The list of addresses of blocks surrounding block 'lb',
!!
!!   logical, intent(in) :: lcc
!!        A logical switch controlling whether 'unk1' or 'work1' data is filled.
!!        I.e. if true then only 'unk' OR 'work' will have guardcells filled 
!!        depending on the value of 'iopt'.
!!  
!!   logical, intent(in) :: lfc
!!        A logical switch controlling whether 'facevarx1, y1, and z' have 
!!        their guardcells filled.
!!
!!   logical, intent(in) :: lec
!!        A logical switch controlling whether 'unk_e_x1, y1, and z1' have 
!!        their guardcells filled.
!!
!!   logical, intent(in) :: lnc
!!        A logical switch controlling whether 'unk_n' has its guardcells filled.
!!
!!   integer, intent(in) :: icoord
!!        An integer switch used to select which faces of the block are to be considered.
!!        If icoord = 0 all faces are considered.  If icoord = 1 only faces perpendicular
!!        to the x-axis are considered. If icoord = 2 only faces perpendicular to the
!!        y-axis are considered.  If icoord = 3 only faces perpendicular to the
!!        z-axis are considered.
!!
!!   logical, intent(in) :: ldiag
!!        A logical switch which controls whether guardcells corresponding to 
!!        neighbor blocks diagonally opposite block edges and corners are filled.
!!
!!   integer, intent(in) :: nlayersx, nlayersy, nlayersz
!!        Arguments to specify the number of guardcells to fill in each 
!!        coordninate direction.
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!   mpif.h
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   timings
!!   workspace
!!   prolong_arrays
!!   paramesh_interfaces
!!
!! CALLS
!!
!!   amr_1blk_cc_prol_gen_work_fun
!!   amr_1blk_cc_prol_gen_unk_fun
!!   amr_1blk_fc_prol_gen_fun
!!   amr_1blk_ec_prol_gen_fun
!!   amr_1blk_nc_prol_gen_fun
!!   amr_1blk_fc_prol_dbz
!!   amr_1blk_fc_clean_divb
!!
!! RETURNS
!!
!!   Does not return anything.  Upon exit, guarcells of the block being operated on 
!!   (which is at a jump in refinement) has its guardcells filled.
!!
!! DESCRIPTION
!!
!!   This routine manages the exchange of guard cell information between
!!   blocks required to fill guard cells on block lb, assuming that 
!!   exchange is only required from surrounding blocks at a coarser 
!!   refinement level.
!!
!!   If you are using an odd number of grid cells, then the interface
!!   condition implicit in this routine is that the finer block uses
!!   data prolonged from its parent on the faces which border coarser
!!   grid blocks.
!!
!! AUTHORS
!!
!!   Peter MacNeice (July 1998) with modifications by 
!!   Kevin Olson (2003) for layered guardcell filling.
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

!#define MLRMOD
#undef MLRMOD

      subroutine mpi_amr_1blk_guardcell_c_to_f( & 
     &                    mype,lb,pe,iopt,nlayers, & 
     &                    surrblks,lcc,lfc,lec,lnc,icoord,ldiag, & 
     &                    nlayersx,nlayersy,nlayersz,ipolar)


      use paramesh_dimensions
      use physicaldata
      use tree
      use timings
      use workspace

      use prolong_arrays, only :  & 
     &     prol_fc_dbz, & 
     &     prol_fc_dbz_ivar, & 
     &     prol_fc_dbz_n,  & 
     &     prol_fc_clean_divb

      use paramesh_interfaces, only : amr_1blk_cc_prol_gen_work_fun, & 
     &                                amr_1blk_cc_prol_gen_unk_fun, & 
     &                                amr_1blk_fc_prol_gen_fun, & 
     &                                amr_1blk_ec_prol_gen_fun, & 
     &                                amr_1blk_nc_prol_gen_fun & 
     &                               ,amr_1blk_fc_prol_dbz, & 
     &                                amr_1blk_fc_clean_divb

      implicit none
      include 'mpif.h'

      integer, intent(in) :: mype,iopt,nlayers,lb,pe,icoord
      integer, intent(in) :: surrblks(:,:,:,:)
      logical, intent(in) :: lcc,lfc,lec,lnc,ldiag
      integer, intent(in) :: nlayersx,nlayersy,nlayersz
      integer, intent(in) :: ipolar(2)


!------------------------------------
! local arrays and variables


      integer,save :: cparent(2)
      integer,save :: awhich_child(1)

      integer,save :: cwhich_child

      logical :: lcoarse_neigh

      integer      :: parent_blk,parent_pe
      integer,save :: remote_blk,remote_pe
      integer      :: i, j, k, ii, jj, kk
      integer      :: igc_off_x, igc_off_y, igc_off_z
      integer      :: ioff, joff, koff
      integer      :: nguard0, nlayers0
      integer      :: jfl, jfu, jface
      integer      :: ia, ib, ja, jb, ka, kb
      integer      :: ionea, ioneb, jonea, joneb, konea, koneb
      integer      :: ia1, ib1, ja1, jb1, ka1, kb1
      integer      :: imod, jmod, kmod, idest
      integer      :: iv1, iv2, iv3, iprol

      double precision :: time1

#ifdef MLRMOD
!--------------------------------------------------

! mlr

      logical :: lbalpro

!--------------------------------------------------

      lbalpro = .true. ! .false.

!--------------------------------------------------
! MLRMOD
#endif 

      if (timing_mpi) then
         time1 = mpi_wtime()
      endif
!------------------------------------

      igc_off_x = gc_off_x
      igc_off_y = gc_off_y
      igc_off_z = gc_off_z

!
! Does this block have any neighbors at lower refinement level?
      lcoarse_neigh = .false.
      do k = 2-k3d,2+k3d
      do j = 2-k2d,2+k2d
      do i = 1,3
         if(surrblks(1,i,j,k).lt.1.and.surrblks(1,i,j,k).gt.-20) & 
     &        lcoarse_neigh = .true.
      enddo
      enddo
      enddo
      if(.not.lcoarse_neigh) return

!------------------------------------

!
! Does current block have a parent?
      if(lnew_parent) then
        cparent(:) = parent(:,lb)
      endif
      
      if(cparent(1).gt.-1) then

!
! Get parent's address
      parent_blk = cparent(1)
      parent_pe  = cparent(2)

      cwhich_child = which_child(lb)


!
! compute the offset in the parent block appropriate for this child
       ioff = mod(cwhich_child-1,2)*nxb/2
       joff = mod((cwhich_child-1)/2,2)*nyb/2
       koff = mod((cwhich_child-1)/4,2)*nzb/2

      endif                               ! end of parent if test

!------------------------------------

      nguard0 = nguard
      nlayers0 = nguard

      if(iopt.ge.2) then 
           nguard0 = nguard_work
           nlayers0 = nlayers
      endif

!
! First deal with block's regular faces

      jfl = 1
      jfu = nfaces
      if(icoord.gt.0) then
        jfl = 1 + 2*(icoord-1)
        jfu = jfl + 1
      endif

! cycle through block faces
      do jface = jfl,jfu


        ia = 1 + nguard0
        ib = nxb + nguard0
        ja = 1 + nguard0*k2d
        jb = nyb + nguard0*k2d
        ka = 1 + nguard0*k3d
        kb = nzb + nguard0*k3d
        ionea = 0
        ioneb = 1
        jonea = 0
        joneb = k2d
        konea = 0
        koneb = k3d

        if(jface.eq.1) then
          remote_blk = surrblks(1,1,2,2)
          remote_pe  = surrblks(2,1,2,2)
          ia   = 1 + nguard0 - nlayersx
          ib   = nguard0 + igc_off_x
        elseif(jface.eq.2) then
          remote_blk = surrblks(1,3,2,2)
          remote_pe  = surrblks(2,3,2,2)
          ia   = 1 + nxb + nguard0 - igc_off_x
          ib   = nxb + 2*nguard0
        elseif(jface.eq.3) then
          remote_blk = surrblks(1,2,1,2)
          remote_pe  = surrblks(2,2,1,2)
          ja   = 1 + nguard0 - nlayersy
          jb   = nguard0 + igc_off_y
        elseif(jface.eq.4) then
          remote_blk = surrblks(1,2,3,2)
          remote_pe  = surrblks(2,2,3,2)
          ja   = 1 + nyb + nguard0 - igc_off_y
          jb   = nyb + 2*nguard0
        elseif(jface.eq.5) then
          remote_blk = surrblks(1,2,2,1)
          remote_pe  = surrblks(2,2,2,1)
          ka   = 1 + nguard0 - nlayersz
          kb   = nguard0 + igc_off_z
        elseif(jface.eq.6) then
          remote_blk = surrblks(1,2,2,3)
          remote_pe  = surrblks(2,2,2,3)
          ka   = 1 + nzb + nguard0 - igc_off_z
          kb   = nzb + 2*nguard0
        endif

!
! If a neighbor exists at this blocks refinement level then fill guardcells
! from its data.

        if(remote_blk.gt.-20.and.remote_blk.lt.0) then


! interpolate(prolongate) data from the parent to the child
          if(iopt.eq.1.and.lcc) then
             call amr_1blk_cc_prol_gen_unk_fun( & 
     &                      unk1(:,:,:,:,2),ia,ib,ja,jb,ka,kb, & 
     &                                  1,ioff,joff,koff,mype, & 
     &                                  lb,parent_pe,parent_blk)
          elseif(iopt.ge.2) then
             call amr_1blk_cc_prol_gen_work_fun(work1(:,:,:,2), & 
     &                                  ia,ib,ja,jb,ka,kb, & 
     &                                  1,ioff,joff,koff,mype, & 
     &                                  lb,parent_pe,parent_blk, & 
     &                                  interp_mask_work(iopt-1))
          endif

          if(lfc.and.iopt.eq.1) then

#ifdef MLRMOD
       if(prol_fc_dbz)then
       do iprol = 1, prol_fc_dbz_n
          iv1 = prol_fc_dbz_ivar(1,iprol)
          iv2 = prol_fc_dbz_ivar(2,iprol)
          iv3 = prol_fc_dbz_ivar(3,iprol)

             call amr_1blk_fc_prol_dbz(  & 
     &            facevarx1(:,:,:,:,2),  & 
     &            facevary1(:,:,:,:,2),  & 
     &            facevarz1(:,:,:,:,2), & 
     &            nfacevar,  & 
     &            iv1,iv2,iv3, & 
     &            ia+ionea,ib+ioneb, & 
     &            ja+jonea,jb+joneb, & 
     &            ka+konea,kb+koneb,1, & 
     &            ioff,joff,koff, & 
     &            mype,lb,parent_pe,parent_blk & 
     &               )
!             call amr_1blk_fc_prol_dbz( 
!     .            recvfx, recvfy, recvfz,nfacevar, 
!     .            iv1,iv2,iv3,
!     .            ia+ionea,ib+ioneb,
!     .            ja+jonea,jb+joneb,
!     .            ka+konea,kb+koneb,1,
!     .            ioff,joff,koff,
!     .            mype,lb,parent_pe,parent_blk
!     .               )

       end do
       end if

!MLRMOD
#endif

            if(jface.eq.2) ionea = 1 
            if(jface.eq.1) ioneb = 0
            call amr_1blk_fc_prol_gen_fun(  & 
     &       facevarx1(:nfacevar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1, & 
     &                 kl_bnd1:ku_bnd1,2), & 
     &                    ia+ionea,ib+ioneb,ja,jb,ka,kb,1, & 
     &                    ioff,joff,koff, & 
     &                    mype,lb,parent_pe,parent_blk,1)

            if(ndim.ge.2) then

            if(jface.eq.4) jonea = 1
            if(jface.eq.3) joneb = 0
            call amr_1blk_fc_prol_gen_fun(  & 
     &       facevary1(:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d, & 
     &                 kl_bnd1:ku_bnd1,2), & 
     &                    ia,ib,ja+jonea,jb+joneb,ka,kb,1, & 
     &                    ioff,joff,koff, & 
     &                    mype,lb,parent_pe,parent_blk,2)
            endif

            if(ndim.eq.3) then

            if(jface.eq.6) konea = 1
            if(jface.eq.5) koneb = 0
            call amr_1blk_fc_prol_gen_fun(  & 
     &       facevarz1(:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1, & 
     &                 kl_bnd1:ku_bnd1+k3d,2), & 
     &                    ia,ib,ja,jb,ka+konea,kb+koneb,1, & 
     &                    ioff,joff,koff, & 
     &                    mype,lb,parent_pe,parent_blk,3)
!+ clean_divb
            if (clean_divb) then
            if(prol_fc_clean_divb)then
               call amr_1blk_fc_clean_divb(  & 
     &              nfacevar, & 
     &              ia,ib, & 
     &              ja,jb, & 
     &              ka,kb,    & 
     &              ionea,ioneb, & 
     &              jonea,joneb, & 
     &              konea,koneb,    & 
     &              idest,ioff,joff,koff, & 
     &              mype,lb,parent_pe,parent_blk & 
     &              )
            end if
            end if ! clean_divb

            endif

          endif

          if (ndim > 1) then
          if(lec.and.iopt.eq.1) then

            ionea = 0
            ioneb = 1
            jonea = 0
            joneb = k2d
            konea = 0
            koneb = k3d
            if(jface.eq.2) ionea = 1
            if(jface.eq.1) ioneb = 0
            if(jface.eq.4) jonea = 1
            if(jface.eq.3) joneb = 0
            if(jface.eq.6) konea = 1
            if(jface.eq.5) koneb = 0


            call amr_1blk_ec_prol_gen_fun(  & 
     &       unk_e_x1(:nvaredge,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d, & 
     &                kl_bnd1:ku_bnd1+k3d,2), & 
     &                    ia,ib,ja+jonea,jb+joneb,ka+konea,kb+koneb,1, & 
     &                    ioff,joff,koff,mype,1)

            call amr_1blk_ec_prol_gen_fun(  & 
     &       unk_e_y1(:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1, & 
     &                kl_bnd1:ku_bnd1+k3d,2), & 
     &                    ia+ionea,ib+ioneb,ja,jb,ka+konea,kb+koneb,1, & 
     &                    ioff,joff,koff,mype,2)

            if (ndim == 3) then

            call amr_1blk_ec_prol_gen_fun(  & 
     &       unk_e_z1(:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d, & 
     &                kl_bnd1:ku_bnd1,2), & 
     &                    ia+ionea,ib+ioneb,ja+jonea,jb+joneb,ka,kb,1, & 
     &                    ioff,joff,koff,mype,3)
            endif ! end if (ndim ==3)

          endif

          end if ! end if (ndim > 1)

          if(lnc.and.iopt.eq.1) then

            ionea = 0
            ioneb = 1
            jonea = 0
            joneb = k2d
            konea = 0
            koneb = k3d

            if(jface.eq.2) ionea = 1
            if(jface.eq.1) ioneb = 0
            if(jface.eq.4) jonea = 1
            if(jface.eq.3) joneb = 0
            if(jface.eq.6) konea = 1
            if(jface.eq.5) koneb = 0
            call amr_1blk_nc_prol_gen_fun(  & 
     &       unk_n1(:nvarcorn,il_bnd1:iu_bnd1+1, & 
     &              jl_bnd1:ju_bnd1+k2d,kl_bnd1:ku_bnd1+k3d,2), & 
     &                    ia+ionea,ib+ioneb,ja+jonea,jb+joneb, & 
     &                    ka+konea,kb+koneb,1, & 
     &                    ioff,joff,koff,mype)

          endif

        endif

      enddo                           ! end of loop over block faces

!------------------------------------

      if(ldiag) then

!------------------------------------
      if(ndim.ge.2) then
      if(icoord.ne.3) then
!
! Now fill from edges along the z axis.

      ia = 1 + nguard0
      ib = nxb + nguard0
      ja = 1 + nguard0*k2d
      jb = nyb + nguard0*k2d
      ka = 1 + nguard0*k3d
      kb = nzb + nguard0*k3d

! Loop over the 4 corners
      do jj = 1,3,2
      do ii = 1,3,2
    
          remote_blk = surrblks(1,ii,jj,2)
          remote_pe  = surrblks(2,ii,jj,2)

          imod = ii/2
          jmod = jj/2*k2d

          ia = 1 + (nguard0 - nlayersx) +  & 
     &             (nxb + nlayersx - igc_off_x )*imod
          ib = nguard0 + igc_off_x +  & 
     &             (nxb + nlayersx - igc_off_x )*imod
          ja = 1 + (nguard0 - nlayersy) +  & 
     &             (nyb + nlayersy - igc_off_y )*jmod
          jb = nguard0 + igc_off_y + & 
     &             (nyb + nlayersy - igc_off_y )*jmod

          ia1 = imod
          ib1 = imod
          ja1 = jmod
          jb1 = jmod


        if(remote_blk.gt.-20.and.remote_blk.lt.0) then

! interpolate(prolongate) data from the parent to the child
          if(iopt.eq.1.and.lcc) then
             call amr_1blk_cc_prol_gen_unk_fun( & 
     &                      unk1(:,:,:,:,2),ia,ib,ja,jb,ka,kb, & 
     &                                    1,ioff,joff,koff,mype, & 
     &                                    lb,parent_pe,parent_blk)
          elseif(iopt.ge.2) then
             call amr_1blk_cc_prol_gen_work_fun(work1(:,:,:,2), & 
     &                                    ia,ib,ja,jb,ka,kb, & 
     &                                    1,ioff,joff,koff,mype, & 
     &                                    lb,parent_pe,parent_blk, & 
     &                                    interp_mask_work(iopt-1))
          endif

          if(lfc.and.iopt.eq.1) then

#ifdef MLRMOD
       if(prol_fc_dbz)then
       do iprol = 1, prol_fc_dbz_n
          iv1 = prol_fc_dbz_ivar(1,iprol)
          iv2 = prol_fc_dbz_ivar(2,iprol)
          iv3 = prol_fc_dbz_ivar(3,iprol)

             call amr_1blk_fc_prol_dbz(  & 
     &            facevarx1(:,:,:,:,2),  & 
     &            facevary1(:,:,:,:,2),  & 
     &            facevarz1(:,:,:,:,2), & 
     &            nfacevar,  & 
     &            iv1,iv2,iv3, & 
     &            ia+ia1,ib+ib1, & 
     &            ja+ja1,jb+jb1, & 
     &            ka,kb+1,1, & 
     &            ioff,joff,koff, & 
     &            mype,lb,parent_pe,parent_blk & 
     &               )

!             call amr_1blk_fc_prol_dbz( 
!     .            recvfx, recvfy, recvfz,nfacevar, 
!     .            iv1,iv2,iv3,
!     .            ia+ia1,ib+ib1,
!     .            ja+ja1,jb+jb1,
!     .            ka,kb+1,1,
!     .            ioff,joff,koff,
!     .            mype,lb,parent_pe,parent_blk
!     .               )

       end do
       end if
!MLRMOD
#endif

            call amr_1blk_fc_prol_gen_fun(  & 
     &       facevarx1(1:nfacevar,il_bnd1:iu_bnd1+1, & 
     &                 jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1,2), & 
     &                    ia+ia1,ib+ib1,ja,jb,ka,kb,1, & 
     &                    ioff,joff,koff, & 
     &                    mype,lb,parent_pe,parent_blk,1)

            if(ndim.ge.2) then

            call amr_1blk_fc_prol_gen_fun(  & 
     &       facevary1(1:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d, & 
     &                 kl_bnd1:ku_bnd1,2), & 
     &                    ia,ib,ja+ja1,jb+jb1,ka,kb,1, & 
     &                    ioff,joff,koff, & 
     &                    mype,lb,parent_pe,parent_blk,2)
            endif

            if(ndim.eq.3) then

            call amr_1blk_fc_prol_gen_fun( & 
     &       facevarz1(1:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1, & 
     &                 kl_bnd1:ku_bnd1+k3d,2), & 
     &                    ia,ib,ja,jb,ka,kb+1,1, & 
     &                    ioff,joff,koff, & 
     &                    mype,lb,parent_pe,parent_blk,3)

!+ clean_divb
       if (clean_divb) then
       if(prol_fc_clean_divb)then
          call amr_1blk_fc_clean_divb(  & 
     &         nfacevar, & 
     &              ia,ib, & 
     &              ja,jb, & 
     &              ka,kb,       & 
     &              ia1,ib1, & 
     &              ja1,jb1, & 
     &              0,1,       & 
     &         idest,ioff,joff,koff, & 
     &         mype,lb,parent_pe,parent_blk & 
     &         )
       end if
       end if  ! clean_divb

            endif

          endif

          if (ndim > 1) then
          if(lec.and.iopt.eq.1) then

            ionea = 0
            ioneb = 1
            jonea = 0
            joneb = k2d
            konea = 0
            koneb = k3d

            if(ii.eq.1) ioneb = 0
            if(ii.eq.3) ionea = 1
            if(jj.eq.1) joneb = 0
            if(jj.eq.3) jonea = 1


            call amr_1blk_ec_prol_gen_fun(  & 
     &       unk_e_x1(:nvaredge,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d, & 
     &                kl_bnd1:ku_bnd1+k3d,2), & 
     &                    ia,ib,ja+jonea,jb+joneb,ka+konea,kb+koneb,1, & 
     &                    ioff,joff,koff,mype,1)

            call amr_1blk_ec_prol_gen_fun(  & 
     &       unk_e_y1(:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1, & 
     &                kl_bnd1:ku_bnd1+k3d,2), & 
     &                    ia+ionea,ib+ioneb,ja,jb,ka+konea,kb+koneb,1, & 
     &                    ioff,joff,koff,mype,2)

            if (ndim == 3) then
            call amr_1blk_ec_prol_gen_fun(  & 
     &       unk_e_z1(:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d, & 
     &                kl_bnd1:ku_bnd1,2), & 
     &                    ia+ionea,ib+ioneb,ja+jonea,jb+joneb,ka,kb,1, & 
     &                    ioff,joff,koff,mype,3)
            end if ! if (ndim == 3)

          endif

          end if ! if (ndim > 1)

          if(lnc.and.iopt.eq.1) then

            call amr_1blk_nc_prol_gen_fun( & 
     &       unk_n1(1:nvarcorn,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d, & 
     &              kl_bnd1:ku_bnd1+k3d,2), & 
     &                    ia+imod,ib+imod,ja+jmod,jb+jmod, & 
     &                    ka,kb+k3d,1, & 
     &                    ioff,joff,koff,mype)

          endif

        endif


      enddo                                   ! end of ii do loop
      enddo                                   ! end of jj do loop


      endif                                   ! end of if test icoord=3
      endif                                   ! end of if test for 2D corners


      if(ndim.eq.3) then
      if(icoord.ne.2) then

!
! Now fill from edges along the y axis.

      ia = 1 + nguard0
      ib = nxb + nguard0
      ja = 1 + nguard0*k2d
      jb = nyb + nguard0*k2d
      ka = 1 + nguard0*k3d
      kb = nzb + nguard0*k3d

! Loop over the 4 corners
      do kk = 1,3,2
      do ii = 1,3,2
    
          remote_blk = surrblks(1,ii,2,kk)
          remote_pe  = surrblks(2,ii,2,kk)

          imod = ii/2
          kmod = kk/2*k3d

          ia = 1 + (nguard0 - nlayersx) +  & 
     &             (nxb + nlayersx - igc_off_x )*imod
          ib = nguard0 + igc_off_x +  & 
     &             (nxb + nlayersx - igc_off_x )*imod
          ka = 1 + (nguard0 - nlayersz) +  & 
     &             (nzb + nlayersz - igc_off_z )*kmod
          kb = nguard0 + igc_off_z + & 
     &             (nzb + nlayersz - igc_off_z )*kmod

          ia1 = imod
          ib1 = imod
          ka1 = kmod
          kb1 = kmod


        if(remote_blk.gt.-20.and.remote_blk.lt.0) then

! interpolate(prolongate) data from the parent to the child
          if(iopt.eq.1.and.lcc) then
             call amr_1blk_cc_prol_gen_unk_fun( & 
     &                      unk1(:,:,:,:,2),ia,ib,ja,jb,ka,kb, & 
     &                                   1,ioff,joff,koff,mype, & 
     &                                   lb,parent_pe,parent_blk)
          elseif(iopt.ge.2) then
             call amr_1blk_cc_prol_gen_work_fun(work1(:,:,:,2), & 
     &                                   ia,ib,ja,jb,ka,kb, & 
     &                                   1,ioff,joff,koff,mype, & 
     &                                   lb,parent_pe,parent_blk, & 
     &                                   interp_mask_work(iopt-1))
          endif

          if(lfc.and.iopt.eq.1) then

#ifdef MLRMOD
       if(prol_fc_dbz)then
       do iprol = 1, prol_fc_dbz_n
          iv1 = prol_fc_dbz_ivar(1,iprol)
          iv2 = prol_fc_dbz_ivar(2,iprol)
          iv3 = prol_fc_dbz_ivar(3,iprol)

             call amr_1blk_fc_prol_dbz(  & 
     &            facevarx1(:,:,:,:,2),  & 
     &            facevary1(:,:,:,:,2),  & 
     &            facevarz1(:,:,:,:,2), & 
     &            nfacevar,  & 
     &            iv1,iv2,iv3, & 
     &            ia+ia1,ib+ib1, & 
     &            ja,jb+1, & 
     &            ka+ka1,kb+kb1,1, & 
     &            ioff,joff,koff, & 
     &            mype,lb,parent_pe,parent_blk & 
     &               )
!             call amr_1blk_fc_prol_dbz( 
!     .            recvfx, recvfy, recvfz,nfacevar, 
!     .            iv1,iv2,iv3,
!     .            ia+ia1,ib+ib1,
!     .            ja,jb+1,
!     .            ka+ka1,kb+kb1,1,
!     .            ioff,joff,koff,
!     .            mype,lb,parent_pe,parent_blk
!     .               )
      end do
      end if
!MLRMOD
#endif

            call amr_1blk_fc_prol_gen_fun(  & 
     &       facevarx1(1:nfacevar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1, & 
     &                 kl_bnd1:ku_bnd1,2), & 
     &                    ia+ia1,ib+ib1,ja,jb,ka,kb,1, & 
     &                    ioff,joff,koff, & 
     &                    mype,lb,parent_pe,parent_blk,1)

            if(ndim.ge.2) then

            call amr_1blk_fc_prol_gen_fun(  & 
     &       facevary1(1:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d, & 
     &                 kl_bnd1:ku_bnd1,2), & 
     &                    ia,ib,ja,jb+1,ka,kb,1, & 
     &                    ioff,joff,koff, & 
     &                    mype,lb,parent_pe,parent_blk,2)
            endif

            if(ndim.eq.3) then

            call amr_1blk_fc_prol_gen_fun(  & 
     &       facevarz1(1:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1, & 
     &                 kl_bnd1:ku_bnd1+k3d,2), & 
     &                    ia,ib,ja,jb,ka+ka1,kb+kb1,1, & 
     &                    ioff,joff,koff, & 
     &                    mype,lb,parent_pe,parent_blk,3)
!+ clean_divb
       if (clean_divb) then
       if(prol_fc_clean_divb)then
          call amr_1blk_fc_clean_divb(  & 
     &         nfacevar, & 
     &              ia,ib, & 
     &              ja,jb, & 
     &              ka,kb,       & 
     &              ia1,ib1, & 
     &              0,1, & 
     &              ka1,kb1, & 
     &         idest,ioff,joff,koff, & 
     &         mype,lb,parent_pe,parent_blk & 
     &         )
       end if
       end if ! clean_divb
            endif

          endif



          if (ndim > 1) then
          if(lec.and.iopt.eq.1) then

            ionea = 0
            ioneb = 1
            jonea = 0
            joneb = k2d
            konea = 0
            koneb = k3d

            if(ii.eq.1) ioneb = 0
            if(ii.eq.3) ionea = 1
            if(kk.eq.1) koneb = 0
            if(kk.eq.3) konea = 1

            call amr_1blk_ec_prol_gen_fun(  & 
     &       unk_e_x1(1:nvaredge,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d, & 
     &                kl_bnd1:ku_bnd1+k3d,2), & 
     &                    ia,ib,ja+jonea,jb+joneb,ka+konea,kb+koneb,1, & 
     &                    ioff,joff,koff,mype,1)

            call amr_1blk_ec_prol_gen_fun(  & 
     &       unk_e_y1(1:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1, & 
     &                kl_bnd1:ku_bnd1+k3d,2), & 
     &                    ia+ionea,ib+ioneb,ja,jb,ka+konea,kb+koneb,1, & 
     &                    ioff,joff,koff,mype,2)

            if (ndim == 3) then

            call amr_1blk_ec_prol_gen_fun(  & 
     &       unk_e_z1(1:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d, & 
     &                kl_bnd1:ku_bnd1,2), & 
     &                    ia+ionea,ib+ioneb,ja+jonea,jb+joneb,ka,kb,1, & 
     &                    ioff,joff,koff,mype,3)
            end if ! end if (ndim == 3)

          endif

          end if ! end  if (ndim > 1)

          if(lnc.and.iopt.eq.1) then

            call amr_1blk_nc_prol_gen_fun(  & 
     &       unk_n1(1:nvarcorn,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d, & 
     &              kl_bnd1:ku_bnd1+k3d,2), & 
     &                    ia+imod,ib+imod,ja,jb+k2d, & 
     &                    ka+kmod,kb+kmod,1, & 
     &                    ioff,joff,koff,mype)

          endif

        endif


      enddo
      enddo

      endif                                    ! end of icoord=2 if test
!
! Now fill from edges along the x axis.

      if(icoord.ne.1) then

      ia = 1 + nguard0
      ib = nxb + nguard0
      ja = 1 + nguard0*k2d
      jb = nyb + nguard0*k2d
      ka = 1 + nguard0*k3d
      kb = nzb + nguard0*k3d

! Loop over the 4 corners
      do kk = 1,3,2
      do jj = 1,3,2
    
          remote_blk = surrblks(1,2,jj,kk)
          remote_pe  = surrblks(2,2,jj,kk)

          jmod = jj/2*k2d
          kmod = kk/2*k3d

          ja = 1 + (nguard0 - nlayersy) +  & 
     &             (nyb + nlayersy - igc_off_y )*jmod
          jb = nguard0 + igc_off_y + & 
     &             (nyb + nlayersy - igc_off_y )*jmod
          ka = 1 + (nguard0 - nlayersz) +  & 
     &             (nzb + nlayersz - igc_off_z )*kmod
          kb = nguard0 + igc_off_z + & 
     &             (nzb + nlayersz - igc_off_z )*kmod

          ja1 = jmod
          jb1 = jmod
          ka1 = kmod
          kb1 = kmod

        if(remote_blk.gt.-20.and.remote_blk.lt.0) then

! interpolate(prolongate) data from the parent to the child
          if(iopt.eq.1.and.lcc) then
             call amr_1blk_cc_prol_gen_unk_fun( & 
     &                      unk1(:,:,:,:,2),ia,ib,ja,jb,ka,kb, & 
     &                                    1,ioff,joff,koff,mype, & 
     &                                    lb,parent_pe,parent_blk)
          elseif(iopt.ge.2) then
             call amr_1blk_cc_prol_gen_work_fun(work1(:,:,:,2), & 
     &                                     ia,ib,ja,jb,ka,kb, & 
     &                                     1,ioff,joff,koff,mype, & 
     &                                     lb,parent_pe,parent_blk, & 
     &                                     interp_mask_work(iopt-1))
          endif

          if(lfc.and.iopt.eq.1) then

#ifdef MLRMOD
       if(prol_fc_dbz)then
       do iprol = 1, prol_fc_dbz_n
          iv1 = prol_fc_dbz_ivar(1,iprol)
          iv2 = prol_fc_dbz_ivar(2,iprol)
          iv3 = prol_fc_dbz_ivar(3,iprol)

             call amr_1blk_fc_prol_dbz(  & 
     &            facevarx1(:,:,:,:,2),  & 
     &            facevary1(:,:,:,:,2),  & 
     &            facevarz1(:,:,:,:,2), & 
     &            nfacevar,  & 
     &            iv1,iv2,iv3, & 
     &            ia,ib+1, & 
     &            ja+ja1,jb+jb1, & 
     &            ka,kb,1, & 
     &            ioff,joff,koff, & 
     &            mype,lb,parent_pe,parent_blk & 
     &               )
!             call amr_1blk_fc_prol_dbz( 
!     .            recvfx, recvfy, recvfz,nfacevar, 
!     .            iv1,iv2,iv3,
!     .            ia,ib+1,
!     .            ja+ja1,jb+jb1,
!     .            ka,kb,1,
!     .            ioff,joff,koff,
!     .            mype,lb,parent_pe,parent_blk
!     .               )
       end do
       end if
!MLRMOD
#endif

            call amr_1blk_fc_prol_gen_fun(  & 
     &       facevarx1(1:nfacevar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1, & 
     &                 kl_bnd1:ku_bnd1,2), & 
     &                    ia,ib+1,ja,jb,ka,kb,1, & 
     &                    ioff,joff,koff, & 
     &                    mype,lb,parent_pe,parent_blk,1)

            if(ndim.ge.2) then
            call amr_1blk_fc_prol_gen_fun(  & 
     &        facevary1(1:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d, & 
     &                  kl_bnd1:ku_bnd1,2), & 
     &                    ia,ib,ja+ja1,jb+jb1,ka,kb,1, & 
     &                    ioff,joff,koff, & 
     &                    mype,lb,parent_pe,parent_blk,2)
            endif

            if(ndim.eq.3) then
            call amr_1blk_fc_prol_gen_fun(  & 
     &       facevarz1(1:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1, & 
     &                 kl_bnd1:ku_bnd1+k3d,2), & 
     &                    ia,ib,ja,jb,ka+ka1,kb+kb1,1, & 
     &                    ioff,joff,koff, & 
     &                    mype,lb,parent_pe,parent_blk,3)
!+ clean_divb
       if (clean_divb) then
       if(prol_fc_clean_divb)then
          call amr_1blk_fc_clean_divb(  & 
     &         nfacevar, & 
     &              ia,ib, & 
     &              ja,jb, & 
     &              ka,kb,       & 
     &              0,1, & 
     &              ja1,jb1, & 
     &              ka1,kb1, & 
     &         idest,ioff,joff,koff, & 
     &         mype,lb,parent_pe,parent_blk & 
     &         )
       end if
       end if ! clean_divb
            endif

          endif

          if (ndim > 1) then
          if(lec.and.iopt.eq.1) then

            ionea = 0
            ioneb = 1
            jonea = 0
            joneb = k2d
            konea = 0
            koneb = k3d

            if(jj.eq.1) joneb = 0
            if(jj.eq.3) jonea = 1
            if(kk.eq.1) koneb = 0
            if(kk.eq.3) konea = 1

            call amr_1blk_ec_prol_gen_fun(  & 
     &       unk_e_x1(1:nvaredge,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d, & 
     &                kl_bnd1:ku_bnd1+k3d,2), & 
     &                    ia,ib,ja+jonea,jb+joneb,ka+konea,kb+koneb,1, & 
     &                    ioff,joff,koff,mype,1)

            call amr_1blk_ec_prol_gen_fun(  & 
     &       unk_e_y1(1:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1, & 
     &                kl_bnd1:ku_bnd1+k3d,2), & 
     &                    ia+ionea,ib+ioneb,ja,jb,ka+konea,kb+koneb,1, & 
     &                    ioff,joff,koff,mype,2)

            if (ndim == 3) then
            call amr_1blk_ec_prol_gen_fun(  & 
     &       unk_e_z1(1:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d, & 
     &                kl_bnd1:ku_bnd1,2), & 
     &                    ia+ionea,ib+ioneb,ja+jonea,jb+joneb,ka,kb,1, & 
     &                    ioff,joff,koff,mype,3)
            end if ! if (ndim == 3)

          endif

          end if ! if (ndim > 1)

          if(lnc.and.iopt.eq.1) then

            call amr_1blk_nc_prol_gen_fun(  & 
     &       unk_n1(1:nvarcorn,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d, & 
     &              kl_bnd1:ku_bnd1+k3d,2), & 
     &                    ia,ib+1,ja+jmod,jb+jmod, & 
     &                    ka+kmod,kb+kmod,1, & 
     &                    ioff,joff,koff,mype)

          endif

        endif


      enddo
      enddo

      endif                                 ! end of if test for icoord=1

      endif                                 ! end of if test for x and y edges

!------------------------------------

! Finally fill corners in 3D.

!------------------------------------
      if(ndim.eq.3) then


! Loop over the 4 corners
      do kk = 1,3,2
      do jj = 1,3,2
      do ii = 1,3,2
    
          remote_blk = surrblks(1,ii,jj,kk)
          remote_pe  = surrblks(2,ii,jj,kk)


          imod = ii/2
          jmod = jj/2
          kmod = kk/2

          ia = 1 + (nguard0 - nlayersx) +  & 
     &             (nxb + nlayersx - igc_off_x )*imod
          ib = nguard0 + igc_off_x +  & 
     &             (nxb + nlayersx - igc_off_x )*imod
          ja = 1 + (nguard0 - nlayersy) +  & 
     &             (nyb + nlayersy - igc_off_y )*jmod
          jb = nguard0 + igc_off_y + & 
     &             (nyb + nlayersy - igc_off_y )*jmod
          ka = 1 + (nguard0 - nlayersz) +  & 
     &             (nzb + nlayersz - igc_off_z )*kmod
          kb = nguard0 + igc_off_z + & 
     &             (nzb + nlayersz - igc_off_z )*kmod

          ia1 = imod
          ib1 = imod
          ja1 = jmod
          jb1 = jmod
          ka1 = kmod
          kb1 = kmod

        if(remote_blk.gt.-20.and.remote_blk.lt.0) then

! interpolate(prolongate) data from the parent to the child
          if(iopt.eq.1.and.lcc) then
             call amr_1blk_cc_prol_gen_unk_fun( & 
     &                      unk1(:,:,:,:,2),ia,ib,ja,jb,ka,kb, & 
     &                                    1,ioff,joff,koff,mype, & 
     &                                    lb,parent_pe,parent_blk)

          elseif(iopt.ge.2) then
             call amr_1blk_cc_prol_gen_work_fun(work1(:,:,:,2), & 
     &                                    ia,ib,ja,jb,ka,kb, & 
     &                                    1,ioff,joff,koff,mype, & 
     &                                    lb,parent_pe,parent_blk, & 
     &                                    interp_mask_work(iopt-1))
          endif

          if(lfc.and.iopt.eq.1) then

#ifdef MLRMOD
       if(prol_fc_dbz)then
       do iprol = 1, prol_fc_dbz_n
          iv1 = prol_fc_dbz_ivar(1,iprol)
          iv2 = prol_fc_dbz_ivar(2,iprol)
          iv3 = prol_fc_dbz_ivar(3,iprol)

                call amr_1blk_fc_prol_dbz(  & 
     &               facevarx1(:,:,:,:,2),  & 
     &               facevary1(:,:,:,:,2),  & 
     &               facevarz1(:,:,:,:,2), & 
     &               nfacevar,     & 
     &               iv1,iv2,iv3, & 
     &               ia+ia1,ib+ib1,ja+ja1,jb+jb1,ka+ka1,kb+kb1,1, & 
     &               ioff,joff,koff, & 
     &               mype,lb,parent_pe,parent_blk    & 
     &               )

!                call amr_1blk_fc_prol_dbz( 
!     .               recvfx, recvfy, recvfz,nfacevar,    
!     .               iv1,iv2,iv3,
!     .               ia+ia1,ib+ib1,ja+ja1,jb+jb1,ka+ka1,kb+kb1,1,
!     .               ioff,joff,koff,
!     .               mype,lb,parent_pe,parent_blk   
!     .               )

       end do
       end if
!MLRMOD
#endif


            call amr_1blk_fc_prol_gen_fun(  & 
     &       facevarx1(:nfacevar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1, & 
     &                 kl_bnd1:ku_bnd1,2), & 
     &                    ia+ia1,ib+ib1,ja,jb,ka,kb,1, & 
     &                    ioff,joff,koff, & 
     &                    mype,lb,parent_pe,parent_blk,1)

            if(ndim.ge.2) then
            call amr_1blk_fc_prol_gen_fun(  & 
     &       facevary1(:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d, & 
     &                 kl_bnd1:ku_bnd1,2), & 
     &                    ia,ib,ja+ja1,jb+jb1,ka,kb,1, & 
     &                    ioff,joff,koff, & 
     &                    mype,lb,parent_pe,parent_blk,2)
            endif

            if(ndim.eq.3) then
            call amr_1blk_fc_prol_gen_fun(  & 
     &       facevarz1(:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1, & 
     &                 kl_bnd1:ku_bnd1+k3d,2), & 
     &                    ia,ib,ja,jb,ka+ka1,kb+kb1,1, & 
     &                    ioff,joff,koff, & 
     &                    mype,lb,parent_pe,parent_blk,3)
!+ clean_divb
        if (clean_divb) then
        if(prol_fc_clean_divb)then
           call amr_1blk_fc_clean_divb(  & 
     &              nfacevar, & 
     &              ia,ib, & 
     &              ja,jb, & 
     &              ka,kb,      & 
     &              ia1,ib1, & 
     &              ja1,jb1, & 
     &              ka1,kb1, & 
     &         idest,ioff,joff,koff, & 
     &         mype,lb,parent_pe,parent_blk & 
     &         )
        end if
        end if ! clean_divb
            endif

          endif


          if (ndim > 1) then
          if(lec.and.iopt.eq.1) then

            ionea = 0
            ioneb = 1
            jonea = 0
            joneb = k2d
            konea = 0
            koneb = k3d

            if(ii.eq.1) ioneb = 0
            if(ii.eq.3) ionea = 1
            if(jj.eq.1) joneb = 0
            if(jj.eq.3) jonea = 1
            if(kk.eq.1) koneb = 0
            if(kk.eq.3) konea = 1

            call amr_1blk_ec_prol_gen_fun(  & 
     &       unk_e_x1(:nvaredge,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d, & 
     &                kl_bnd1:ku_bnd1+k3d,2), & 
     &                    ia,ib,ja+jonea,jb+joneb,ka+konea,kb+koneb,1, & 
     &                    ioff,joff,koff,mype,1)

            call amr_1blk_ec_prol_gen_fun(  & 
     &       unk_e_y1(:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1, & 
     &                kl_bnd1:ku_bnd1+k3d,2), & 
     &                    ia+ionea,ib+ioneb,ja,jb,ka+konea,kb+koneb,1, & 
     &                    ioff,joff,koff,mype,2)

            if (ndim == 3) then
            call amr_1blk_ec_prol_gen_fun(  & 
     &       unk_e_z1(:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d, & 
     &                kl_bnd1:ku_bnd1,2), & 
     &                    ia+ionea,ib+ioneb,ja+jonea,jb+joneb,ka,kb,1, & 
     &                    ioff,joff,koff,mype,3)
            end if ! if (ndim == 3)

          endif

          end if ! if (ndim > 1)

          if(lnc.and.iopt.eq.1) then

            call amr_1blk_nc_prol_gen_fun(  & 
     &       unk_n1(:nvarcorn,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d, & 
     &              kl_bnd1:ku_bnd1+k3d,2), & 
     &                    ia+imod,ib+imod,ja+jmod,jb+jmod, & 
     &                    ka+kmod,kb+kmod,1, & 
     &                    ioff,joff,koff,mype)

          endif

        endif


      enddo
      enddo
      enddo

      endif                                 ! end of loop over 3D corners

!------------------------------------

      endif                                 ! end of ldiag if test

!------------------------------------
      if (timing_mpi) then
       timer_amr_1blk_guardcell_c_to_f =   & 
     &     timer_amr_1blk_guardcell_c_to_f + mpi_wtime() - time1
      endif

      return
      end subroutine mpi_amr_1blk_guardcell_c_to_f


