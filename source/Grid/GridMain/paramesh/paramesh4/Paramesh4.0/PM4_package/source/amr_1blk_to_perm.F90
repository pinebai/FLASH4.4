!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"


      subroutine amr_1blk_to_perm( lcc,lfc,lec,lnc,lb,iopt,idest)




!------------------------------------------------------------------------
!
! This routine copies data from the 1-block working arrays with guardcells
! to the permanent data arrays, which may or may not have permanent
! guardcells, depending on whether NO_PERMANENT_GUARDCELLS is defined 
! in physicaldata.fh.
!
!
! Arguments :
!      lcc          logical       copies cell centered data if true
!      lfc          logical       copies cell face-centered data if true
!      lec          logical       copies cell edge-centered data if true
!      lnc          logical       copies cell corner data if true
!      lb           integer       block into which data is to be copied
!      iopt         integer       data structure to be copied
!      idest        integer       sets value for last dimension index
!                                  in the 1-blk data arrays
!
!
! Written :     Peter MacNeice          February 1999
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use timings
      use workspace

      implicit none

!------------------------------------

      integer, intent(in) :: lb,iopt,idest
      logical, intent(in) :: lcc,lfc,lec,lnc

      integer :: iopt0, ivar, ivar_next
      integer :: nguard0, nguard_work0

      include 'mpif.h'
      double precision :: time1

!------------------------------------


      if (timing_mpi) then
         time1 = mpi_wtime()
      endif

      nguard0 = nguard*(1-npgs)
      nguard_work0 = nguard_work*(1-npgs)


! cell-centered data
       if(lcc) then

         if(iopt.eq.1) then

           if(ngcell_on_cc < nvar) then
             do ivar=1,ngcell_on_cc
               ivar_next = gcell_on_cc_pointer(ivar)
               unk(ivar_next,il_bnd:iu_bnd,jl_bnd:ju_bnd, & 
     &                                     kl_bnd:ku_bnd,lb) & 
     &         = unk1(ivar_next,il_bnd+nguard0:iu_bnd+nguard0, & 
     &               jl_bnd+nguard0*k2d:ju_bnd+nguard0*k2d, & 
     &               kl_bnd+nguard0*k3d:ku_bnd+nguard0*k3d,idest)
             enddo
           else
             unk(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,lb) & 
     &         = unk1(:,il_bnd+nguard0:iu_bnd+nguard0, & 
     &               jl_bnd+nguard0*k2d:ju_bnd+nguard0*k2d, & 
     &               kl_bnd+nguard0*k3d:ku_bnd+nguard0*k3d,idest)
           endif

         elseif(iopt.ge.2) then
          iopt0 = iopt-1
           work(ilw:iuw,jlw:juw,klw:kuw,lb,iopt0) & 
     &       = work1(ilw+nguard_work0:iuw+nguard_work0, & 
     &               jlw+nguard_work0*k2d:juw+nguard_work0*k2d, & 
     &               klw+nguard_work0*k3d:kuw+nguard_work0*k3d, & 
     &               idest)
         endif                           ! end of iopt if test
       endif                             ! end of lcc if test



! cell face-centered data
       if(lfc) then

! x-face

         if(ngcell_on_fc(1) < nfacevar) then
           do ivar=1,ngcell_on_fc(1)
             ivar_next = gcell_on_fc_pointer(1,ivar)
             facevarx(ivar_next,il_bnd:iu_bnd+1, & 
     &                jl_bnd:ju_bnd,kl_bnd:ku_bnd,lb) & 
     &         = facevarx1(ivar_next,il_bnd+nguard0:iu_bnd+nguard0+1, & 
     &                jl_bnd+nguard0*k2d:ju_bnd+nguard0*k2d, & 
     &                kl_bnd+nguard0*k3d:ku_bnd+nguard0*k3d,idest)
           enddo
         else
           facevarx(1:nfacevar,il_bnd:iu_bnd+1, & 
     &                       jl_bnd:ju_bnd,kl_bnd:ku_bnd,lb) & 
     &       = facevarx1(1:nfacevar,il_bnd+nguard0:iu_bnd+nguard0+1, & 
     &               jl_bnd+nguard0*k2d:ju_bnd+nguard0*k2d, & 
     &               kl_bnd+nguard0*k3d:ku_bnd+nguard0*k3d,idest)
         endif

         if (ndim > 1) then
! y-face
         if(ngcell_on_fc(2) < nfacevar) then
           do ivar=1,ngcell_on_fc(2)
             ivar_next = gcell_on_fc_pointer(2,ivar)
             facevary(ivar_next,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d, & 
     &                         kl_bnd:ku_bnd,lb) & 
     &       = facevary1(ivar_next,il_bnd+nguard0:iu_bnd+nguard0, & 
     &               jl_bnd+nguard0*k2d:ju_bnd+nguard0*k2d+k2d, & 
     &               kl_bnd+nguard0*k3d:ku_bnd+nguard0*k3d,idest)
           enddo
         else
           facevary(1:nfacevar,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d, & 
     &                         kl_bnd:ku_bnd,lb) & 
     &       = facevary1(1:nfacevar,il_bnd+nguard0:iu_bnd+nguard0, & 
     &               jl_bnd+nguard0*k2d:ju_bnd+nguard0*k2d+k2d, & 
     &               kl_bnd+nguard0*k3d:ku_bnd+nguard0*k3d,idest)
         endif

         end if

         if (ndim == 3) then
! z-face
         if(ngcell_on_fc(3) < nfacevar) then
           do ivar=1,ngcell_on_fc(3)
             ivar_next = gcell_on_fc_pointer(3,ivar)
             facevarz(ivar_next,il_bnd:iu_bnd,jl_bnd:ju_bnd, & 
     &                         kl_bnd:ku_bnd+k3d,lb) & 
     &       = facevarz1(ivar_next,il_bnd+nguard0:iu_bnd+nguard0, & 
     &               jl_bnd+nguard0*k2d:ju_bnd+nguard0*k2d, & 
     &               kl_bnd+nguard0*k3d:ku_bnd+nguard0*k3d+k3d,idest)
           enddo
         else
           facevarz(1:nfacevar,il_bnd:iu_bnd,jl_bnd:ju_bnd, & 
     &                         kl_bnd:ku_bnd+k3d,lb) & 
     &       = facevarz1(1:nfacevar,il_bnd+nguard0:iu_bnd+nguard0, & 
     &               jl_bnd+nguard0*k2d:ju_bnd+nguard0*k2d, & 
     &               kl_bnd+nguard0*k3d:ku_bnd+nguard0*k3d+k3d,idest)
         endif

         end if

        endif                     ! end of lfc if test


! cell edge-centered data

       if (ndim > 1) then
       if(lec) then
! x-edge
         if(ngcell_on_ec(1) < nvaredge) then
           do ivar=1,ngcell_on_ec(1)
             ivar_next = gcell_on_ec_pointer(1,ivar)
             unk_e_x(ivar_next,il_bnd:iu_bnd, & 
     &               jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd+k3d,lb) & 
     &       = unk_e_x1(ivar_next,il_bnd+nguard0:iu_bnd+nguard0, & 
     &               jl_bnd+nguard0*k2d:ju_bnd+(nguard0+1)*k2d, & 
     &               kl_bnd+nguard0*k3d:ku_bnd+(nguard0+1)*k3d,idest)
           enddo
         else
           unk_e_x(1:nvaredge,il_bnd:iu_bnd, & 
     &               jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd+k3d,lb) & 
     &       = unk_e_x1(1:nvaredge,il_bnd+nguard0:iu_bnd+nguard0, & 
     &               jl_bnd+nguard0*k2d:ju_bnd+(nguard0+1)*k2d, & 
     &               kl_bnd+nguard0*k3d:ku_bnd+(nguard0+1)*k3d,idest)
         endif
         if(ngcell_on_ec(2) < nvaredge) then
           do ivar=1,ngcell_on_ec(2)
             ivar_next = gcell_on_ec_pointer(2,ivar)
             unk_e_y(ivar_next,il_bnd:iu_bnd+1, & 
     &               jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d,lb) & 
     &       = unk_e_y1(ivar_next,il_bnd+nguard0:iu_bnd+1+nguard0, & 
     &               jl_bnd+nguard0*k2d:ju_bnd+nguard0*k2d, & 
     &               kl_bnd+nguard0*k3d:ku_bnd+(nguard0+1)*k3d,idest)

           enddo
         else
           unk_e_y(1:nvaredge,il_bnd:iu_bnd+1, & 
     &               jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d,lb) & 
     &       = unk_e_y1(1:nvaredge,il_bnd+nguard0:iu_bnd+1+nguard0, & 
     &               jl_bnd+nguard0*k2d:ju_bnd+nguard0*k2d, & 
     &               kl_bnd+nguard0*k3d:ku_bnd+(nguard0+1)*k3d,idest)
         endif

         if (ndim == 3) then
! z-edge
         if(ngcell_on_ec(3) < nvaredge) then
           do ivar=1,ngcell_on_ec(3)
             ivar_next = gcell_on_ec_pointer(3,ivar)
             unk_e_z(ivar_next,il_bnd:iu_bnd+1, & 
     &               jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd,lb) & 
     &       = unk_e_z1(ivar_next,il_bnd+nguard0:iu_bnd+1+nguard0, & 
     &               jl_bnd+nguard0*k2d:ju_bnd+(nguard0+1)*k2d, & 
     &               kl_bnd+nguard0*k3d:ku_bnd+nguard0*k3d,idest)

           enddo
         else
           unk_e_z(1:nvaredge,il_bnd:iu_bnd+1, & 
     &               jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd,lb) & 
     &       = unk_e_z1(1:nvaredge,il_bnd+nguard0:iu_bnd+1+nguard0, & 
     &               jl_bnd+nguard0*k2d:ju_bnd+(nguard0+1)*k2d, & 
     &               kl_bnd+nguard0*k3d:ku_bnd+nguard0*k3d,idest)
         endif

         end if

        endif                     ! end of lec if test

        end if

! cell corner data
       if(lnc) then
         if(ngcell_on_nc < nvarcorn) then
           do ivar=1,ngcell_on_nc
             ivar_next = gcell_on_nc_pointer(ivar)
             unk_n(ivar_next,il_bnd:iu_bnd+1, & 
     &               jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd+k3d,lb) & 
     &       = unk_n1(ivar_next,il_bnd+nguard0:iu_bnd+1+nguard0, & 
     &               jl_bnd+nguard0*k2d:ju_bnd+(nguard0+1)*k2d, & 
     &               kl_bnd+nguard0*k3d:ku_bnd+(nguard0+1)*k3d,idest)
           enddo
         else
           unk_n(1:nvarcorn,il_bnd:iu_bnd+1, & 
     &               jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd+k3d,lb) & 
     &       = unk_n1(1:nvarcorn,il_bnd+nguard0:iu_bnd+1+nguard0, & 
     &               jl_bnd+nguard0*k2d:ju_bnd+(nguard0+1)*k2d, & 
     &               kl_bnd+nguard0*k3d:ku_bnd+(nguard0+1)*k3d,idest)
         endif
        endif                     ! end of lnc if test


        if (timing_mpi) then
              timer_amr_1blk_to_perm(iopt) =   & 
     &                           timer_amr_1blk_to_perm(iopt) & 
     &                          + mpi_wtime() - time1
        endif
      return
      end subroutine amr_1blk_to_perm
