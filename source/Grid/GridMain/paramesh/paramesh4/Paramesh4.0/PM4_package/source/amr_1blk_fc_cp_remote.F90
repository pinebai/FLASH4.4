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


      subroutine amr_1blk_fc_cp_remote(mype,remote_pe,remote_block, & 
     &   idest,id,jd,kd,is,js,ks,ilays,jlays,klays,ip1,jp1,kp1, & 
     &   ip2,jp2,kp2,iface,nblk_ind,ipolar)


!------------------------------------------------------------------------
!
! This routine copies guard cell information for cell face centered
! data to face iface in block
! idest, from the appropriate face of the neighboring block, assuming
! that the neighboring block is on a different processor.
! This routine is similar to amr_face_cp_remote but allows more flexible
! control over the list of cells which is filled.
! It can be easily edited to alter the data pattern required for schemes
! of different order.
!
! Arguments:
!      mype             local processor
!      remote_pe        remote processor
!      remote_block     local block id of the block to be copied from
!                        the remote processor
!      idest            selects the storage space in data_1blk.fh which is to
!                        be used in this call. If the leaf node is having its
!                        guardcells filled then set this to 1, if its parent
!                        is being filled set it to 2.
!      id               lower limit of index range of points in x direction
!                        on destination block
!      jd               lower limit of index range of points in y direction
!                        on destination block
!      kd               lower limit of index range of points in z direction
!                        on destination block
!      is               lower limit of index range of points in x direction
!                        on source block
!      js               lower limit of index range of points in y direction
!                        on source block
!      ks               lower limit of index range of points in z direction
!                        on source block
!      ilay             no. of mesh points in x direction to be copied
!      jlay             no. of mesh points in y direction to be copied
!      klay             no. of mesh points in z direction to be copied
!      ip1              offset added to index range defined by (id,ilay)
!                        0 if guardcells are at lower end of i index
!                        1 if guardcells are at upper end of i index
!      jp1              offset added to index range defined by (jd,jlay)
!                        0 if guardcells are at lower end of j index
!                        1 if guardcells are at upper end of j index
!      kp1              offset added to index range defined by (kd,klay)
!                        0 if guardcells are at lower end of k index
!                        1 if guardcells are at upper end of k index
!      ip2              extend range in i coord for facevarx by this amount
!                        must be set to either 1 or 0
!      jp2              extend range in j coord for facevary by this amount
!                        must be set to either 1 or 0
!      kp2              extend range in k coord for facevarz by this amount
!                        must be set to either 1 or 0
!      iface            contains the block face on input. If this is
!                        set between 1 to 6 the facevar variables on this
!                        face are averaged with those on a neighbor
!                        at the same refinement level. If iface=0 this
!                        averaging is not done.
!
!
!
! Written :     Peter MacNeice          July 1998
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton
      use paramesh_interfaces, only : amr_mpi_find_blk_in_buffer
      use paramesh_mpi_interfaces, only : mpi_set_message_limits

      implicit none

!-------------------------

      integer, intent(in) :: mype,remote_pe,remote_block
      integer, intent(in) :: idest,id,jd,kd,is,js,ks
      integer, intent(in) :: ilays,jlays,klays
      integer, intent(in) :: ip1,jp1,kp1,ip2,jp2,kp2,iface
      integer, intent(in) :: nblk_ind
      integer, intent(in) :: ipolar(:)
!-------------------------
! local variables

      real, allocatable :: recvxf(:,:,:,:),  & 
     &                     recvyf(:,:,:,:),  & 
     &                     recvzf(:,:,:,:)
#ifndef LIBRARY
      save recvxf,recvyf,recvzf
#endif

      integer :: il,jl,kl,id1,jd1,kd1,is1,js1,ks1,js2
      integer :: ilo,ihi,jlo,jhi,klo,khi
      integer :: ill,jll,kll
      integer :: index0, ii, jj, kk, i, j, k, jbface
      integer :: ia, ib, ja, jb, ka, kb
      integer :: ivar, ivar_next
      integer :: nguard0, jstride
      logical, save :: first_call = .true.
      integer :: dtype
      integer :: vtype
      logical :: lfound
!-------------------------

      nguard0 = nguard*npgs


      if (force_consistency) then

#ifndef LIBRARY
      if (first_call) then
#endif

      if (no_permanent_guardcells) then
       allocate(recvxf(nbndvar,il_bnd:iu_bnd+1,jl_bnd:ju_bnd, & 
     &                 kl_bnd:ku_bnd))
       allocate(recvyf(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d, & 
     &                 kl_bnd:ku_bnd))
       allocate(recvzf(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd, & 
     &                 kl_bnd:ku_bnd+k3d))
      else
       allocate(recvxf(nbndvar,1:2,jl_bnd:ju_bnd,kl_bnd:ku_bnd))
       allocate(recvyf(nbndvar,il_bnd:iu_bnd,1:1+k2d,kl_bnd:ku_bnd))
       allocate(recvzf(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,1:1+k3d))
      end if

      recvxf = 0.
      recvyf = 0.
      recvzf = 0.

#ifndef LIBRARY
      first_call = .false.
      endif
#endif
      endif


!
! Adjust index ranges
      il = ilays-1 
      jl = (jlays-1)*k2d 
      kl = (klays-1)*k3d 

      id1 = id + ip1
      jd1 = jd + jp1*k2d
      kd1 = kd + kp1*k3d
      is1 = is + ip1
      js1 = js + jp1*k2d
      ks1 = ks + kp1*k3d
      ja  = js
      jb  = js + jl 
      js2 = js1 + jl + jp2
      jstride = 1

      if(spherical_pm) then
      if(lsingular_line) then
      if(ipolar(1).eq.-1.and.jd.le.nguard) then
        jstride = -1
        if (no_permanent_guardcells) then
          ja = nguard-jd+1
          jb = ja - jl
        else
          ja = 2*nguard
          jb = ja -nguard + k2d
          js1 = 2*(nguard+1) - jd1
          js2 = js1 + jl*jstride
        endif
      elseif(ipolar(2).eq.+1.and.jd.gt.nyb+nguard) then
        jstride = -1
        if (no_permanent_guardcells) then
          ja = (nyb+1)-(jd-(nyb+nguard))
          jb = ja - jl
          js1 = (nyb+1)-(jd1-(nyb+nguard+1))
        else
          ja = (nyb+nguard)-(jd-(nyb+nguard+1))
          jb = (nyb+nguard)-((jd+jl)-(nyb+nguard+1))
          js1 = (nyb+nguard+1)-(jd1-(nyb+nguard+1))
          js2 = js1 + jl*jstride
        endif
      endif
      endif
      endif

!--
      if(remote_block.le.lnblocks.and.remote_pe.eq.mype) then
!--

       if (no_permanent_guardcells) then
       if(.not.l_f_to_c) then
!--
!       write(*,454) id,jd,kd,il,jl,kl
454   format('range : i-kd ',3(1x,i2),' i-kl ',3(1x,i2))
!       write(*,455) jd,jd+jl,ja,jb,jstride,ipolar
455   format('FXNPG : j ',2(1x,i2),' js ',2(1x,i2),' polar ',2(1x,i2))
       facevarx1(1:nbndvar,id1:id1+il+ip2,jd:jd+jl,kd:kd+kl,idest) & 
     &    = gt_facevarx(1:nbndvar,is1:is1+il+ip2, & 
     &                            ja:jb:jstride, & 
     &                            ks:ks+kl, & 
     &                  remote_block)

       if(ndim.ge.2) then
!       write(*,456) jd1,jd1+jl+jp2,js1,js1+(jl+jp2)*jstride,
!     .               jstride,ipolar
456   format('FYNPG : j ',2(1x,i2),' js1 ',3(1x,i2),' polar ',2(1x,i2))
         facevary1(1:nbndvar,id:id+il,jd1:jd1+jl+jp2,kd:kd+kl,idest) & 
     &      = gt_facevary(1:nbndvar,is:is+il, & 
     &                              js1:js1+(jl+jp2)*jstride:jstride, & 
     &                              ks:ks+kl, & 
     &                    remote_block)
       endif

       if(ndim.eq.3) then
         facevarz1(1:nbndvar,id:id+il,jd:jd+jl,kd1:kd1+kl+kp2,idest) & 
     &      = gt_facevarz(1:nbndvar,is:is+il, & 
     &                              ja:jb:jstride, & 
     &                              ks1:ks1+kl+kp2, & 
     &                    remote_block)
       endif
!--
       else
!--
       facevarx1_fl(1:nbndvar,id1:id1+il+ip2,jd:jd+jl,kd:kd+kl) & 
     &    = gt_facevarx(1:nbndvar,is1:is1+il+ip2, & 
     &                            ja:jb:jstride, & 
     &                            ks:ks+kl, & 
     &                  remote_block)

       if(ndim.ge.2) then
         facevary1_fl(1:nbndvar,id:id+il,jd1:jd1+jl+jp2,kd:kd+kl) & 
     &      = gt_facevary(1:nbndvar,is:is+il, & 
     &                              js1:js1+(jl+jp2)*jstride:jstride, & 
     &                              ks:ks+kl, & 
     &                    remote_block)
       endif

       if(ndim.eq.3) then
         facevarz1_fl(1:nbndvar,id:id+il,jd:jd+jl,kd1:kd1+kl+kp2) & 
     &      = gt_facevarz(1:nbndvar,is:is+il, & 
     &                              ja:jb:jstride, & 
     &                              ks1:ks1+kl+kp2, & 
     &                    remote_block)
       endif

       endif

       if (force_consistency) then
#ifndef LIBRARY
#ifdef NO_PERMANENT_GUARDCELLS
      if(.not.l_f_to_c) then
       recvxf(:,1,:,:) = gt_facevarx(:,1,:,:,remote_block)
       recvxf(:,2,:,:) = gt_facevarx(:,nxb+1,:,:,remote_block)
       if(ndim.ge.2) then
       recvyf(:,:,1,:) = gt_facevary(:,:,1,:,remote_block)
       recvyf(:,:,1+k2d,:) = gt_facevary(:,:,nyb+k2d,:,remote_block)
       endif
       if(ndim.eq.3) then
       recvzf(:,:,:,1) = gt_facevarz(:,:,:,1,remote_block)
       recvzf(:,:,:,1+k3d) = gt_facevarz(:,:,:,nzb+k3d,remote_block)
       endif
      endif
#endif  /* NO_PERMANENT_GUARDCELLS */
#else   /* LIBRARY */
      if(.not.l_f_to_c) then
       recvxf(:,1,:,:) = gt_facevarx(:,1,:,:,remote_block)
       recvxf(:,2,:,:) = gt_facevarx(:,nxb+1,:,:,remote_block)
       if(ndim.ge.2) then
       recvyf(:,:,1,:) = gt_facevary(:,:,1,:,remote_block)
       recvyf(:,:,1+k2d,:) = gt_facevary(:,:,nyb+k2d,:,remote_block)
       endif
       if(ndim.eq.3) then
       recvzf(:,:,:,1) = gt_facevarz(:,:,:,1,remote_block)
       recvzf(:,:,:,1+k3d) = gt_facevarz(:,:,:,nzb+k3d,remote_block)
       endif
      endif
#endif  /* LIBRARY */
      endif

      else ! no_permanent_guardcells

       if(.not.l_f_to_c) then
!--
       facevarx1(1:nbndvar,id1:id1+il+ip2,jd:jd+jl,kd:kd+kl,idest) & 
     &    = facevarx(1:nbndvar,is1:is1+il+ip2,ja:jb:jstride, & 
     &                         ks:ks+kl, & 
     &               remote_block)

       if(ndim.ge.2) then
         facevary1(1:nbndvar,id:id+il,jd1:jd1+jl+jp2,kd:kd+kl,idest) & 
     &      = facevary(1:nbndvar,is:is+il, & 
     &                           js1:js2:jstride, & 
     &                           ks:ks+kl, & 
     &                 remote_block)
       endif

       if(ndim.eq.3) then
         facevarz1(1:nbndvar,id:id+il,jd:jd+jl,kd1:kd1+kl+kp2,idest) & 
     &      = facevarz(1:nbndvar,is:is+il,ja:jb:jstride,ks1:ks1+kl+kp2, & 
     &                 remote_block)
       endif
!--
       else
!--
       facevarx1_fl(1:nbndvar,id1:id1+il+ip2,jd:jd+jl,kd:kd+kl) & 
     &    =  facevarx(1:nbndvar,is1:is1+il+ip2,ja:jb:jstride, & 
     &                          ks:ks+kl, & 
     &                remote_block)

       if(ndim.ge.2) then
         facevary1_fl(1:nbndvar,id:id+il,jd1:jd1+jl+jp2,kd:kd+kl) & 
     &      =    facevary(1:nbndvar,is:is+il, & 
     &                              js1:js1+(jl+jp2)*jstride:jstride, & 
     &                              ks:ks+kl, & 
     &                    remote_block)
       endif

       if(ndim.eq.3) then
         facevarz1_fl(1:nbndvar,id:id+il,jd:jd+jl,kd1:kd1+kl+kp2) & 
     &      =    facevarz(1:nbndvar,is:is+il,ja:jb:jstride, & 
     &                                       ks1:ks1+kl+kp2, & 
     &                    remote_block)
       endif

       endif

      if (force_consistency) then
      if(.not.l_f_to_c) then
       recvxf(:,1,:,:) = gt_facevarx(:,1,:,:,remote_block)
       recvxf(:,2,:,:) = gt_facevarx(:,2,:,:,remote_block)
       if(ndim.ge.2) then
       recvyf(:,:,1,:) = gt_facevary(:,:,1,:,remote_block)
       recvyf(:,:,1+k2d,:) = gt_facevary(:,:,1+k2d,:,remote_block)
       endif
       if(ndim.eq.3) then
       recvzf(:,:,:,1) = gt_facevarz(:,:,:,1,remote_block)
       recvzf(:,:,:,1+k3d) = gt_facevarz(:,:,:,1+k3d,remote_block)
       endif
      endif
      endif

      endif ! no_permanent_guardcells

!--
      else                          ! otherwise if block is remote
!--

        call amr_mpi_find_blk_in_buffer(mype,remote_block, & 
     &                 remote_pe,idest,dtype,index0,lfound)


! If this routine is executing a copy to fill guardcells of a
! leaf blocks^s parent, and the remote block is not found, then
! it is assumed that it is not in the list of buffered remote blocks
! because it is not really needed. Therefore in this case we
! return without copying anything.
        if(idest.eq.2.and.(.not.lfound)) return

! starting index if cell-centered data is also included in recv_buf
        if(l_datapacked(2)) index0 =  & 
     &                      index0 + ngcell_on_cc*message_size_cc(dtype)

        if (l_f_to_c) then
           if (ilays == 2*nguard) then
              ill = nguard
           else
              ill = (ilays-1)/2
           end if
           if (jlays == 2*nguard) then
              jll = nguard
           else
              jll = (jlays-1)/2
           end if
           if (klays == 2*nguard) then
              kll = nguard
           else
              kll = (klays-1)/2
           end if
        else
           ill = ilays
           jll = jlays
           kll = klays
        end if

        vtype = 2
        call mpi_set_message_limits( & 
             &               dtype,ia,ib,ja,jb,ka,kb,vtype, & 
             &               ill,jll,kll)

        ilo = is1
        ihi = is1+il+ip2
        if(iface.eq.1) ihi = nxb+1+nguard0
        if(iface.eq.2) ilo = 1+nguard0


        kk = kd
        do k = ka,kb
           jj = jd
           jstride = 1

           js2 = js
           js1 = js+jl

           if(spherical_pm) then
              if(lsingular_line) then
                 if(ipolar(1).eq.-1.and.jd.le.nguard) then
                    jstride = -1
                    if (no_permanent_guardcells) then
                       js2 = jd + jl + (nguard - 2*(jd+jl)) +1
                       js1 = jd      + (nguard - 2* jd    ) +1
                       jj  = jd + jl
                    else
                       js2 = jd + jl + 2*(nguard - (jd+jl)) +1
                       js1 = jd      + 2*(nguard -  jd    ) +1
                       jj  = jd + jl
                    endif
                 elseif(ipolar(2).eq.+1.and.jd.gt.nyb+nguard) then
                    jstride = -1
                    if (no_permanent_guardcells) then
                       js1 = nyb - ( jd    -(nyb+nguard+1))
                       js2 = nyb - ((jd+jl)-(nyb+nguard+1))
                       jj  = jd + jl
                    else
                       js1 = (nyb+nguard)-( jd    -(nyb+nguard+1))
                       js2 = (nyb+nguard)-((jd+jl)-(nyb+nguard+1))
                       jj  = jd + jl
                    endif
                 endif
              endif
           endif

           do j = ja,jb
              ii = id1
              if(iface.eq.2) ii = id1-1
              do i = ia,ib
                 if (k >= ks .and. k <= ks + kl) then
                    if (j >= js2 .and. j <= js1) then
                       if (i >= ilo .and. i <= ihi) then

                          do ivar=1,ngcell_on_fc(1)
                             ivar_next = gcell_on_fc_pointer(1,ivar)

                             if (.not.l_f_to_c) then
                                if(i.eq.1+nguard0.and.iface.eq.2) then
                                   if (force_consistency) then
                                      recvxf(ivar_next,1,j,k) =  & 
                                           &              temprecv_buf(index0+ivar)
                                   endif
                                elseif(i.eq.nxb+1+nguard0.and.iface.eq.1) then
                                   if (force_consistency) then
                                      recvxf(ivar_next,2,j,k) =  & 
                                           &              temprecv_buf(index0+ivar)
                                   endif
                                else
                                   facevarx1(ivar_next,ii,jj,kk,idest) =  & 
                                        &               temprecv_buf(index0+ivar)
                                endif
                             else
                                facevarx1_fl(ivar_next,ii,jj,kk) =  & 
                                     &               temprecv_buf(index0+ivar)
                             endif

                          enddo

                       end if
                    end if
                 end if
                 if (i >= ilo .and. i <= ihi) ii = ii + 1
                 index0 = index0+ngcell_on_fc(1)

              enddo
              if (j >= js2 .and. j <= js1) jj = jj + jstride
           enddo
           if (k >= ks .and. k <= ks + kl) kk = kk + 1
        enddo

        if(ndim.ge.2) then
           vtype = 3 
           call mpi_set_message_limits( & 
                &               dtype,ia,ib,ja,jb,ka,kb,vtype, & 
                &               ill,jll,kll)


           ! reset js1, js2 to values set on entry
           js1 = js + jp1*k2d
           js2 = js1 + jl + jp2

           jlo = js1
           jhi = js1+jl+jp2
           if(iface.eq.3) jhi = nyb+1+nguard0*k2d
           if(iface.eq.4) jlo = 1 + nguard0*k2d

           kk = kd
           do k = ka,kb
              jj = jd1
              if(iface.eq.4) jj = jd1-1
              jstride = 1

              if(spherical_pm) then
                 if(lsingular_line) then
                    if(ipolar(1).eq.-1.and.jd.le.nguard) then
                       jstride = -1
                       if (no_permanent_guardcells) then
                          jlo = jd + jl + (nguard+2 - 2*(jd+jl))
                          jhi = jd      + (nguard+2 - 2* jd    )
                          jj  = jd + jl
                       else
                          jlo = jd + jl + 2*(nguard+1 - (jd+jl))
                          jhi = jd      + 2*(nguard+1 -  jd    ) 
                          jj  = jd + jl
                       endif
                    elseif(ipolar(2).eq.+1.and.jd.gt.nyb+nguard) then
                       jstride = -1
                       if (no_permanent_guardcells) then
                          jhi = nyb - ( jd    -(nyb+nguard+1))
                          jlo = nyb - ((jd+jl)-(nyb+nguard+1))
                          jj  = jd + jl + 1
                       else
                          jhi = (nyb+nguard)-( jd    -(nyb+nguard+1))
                          jlo = (nyb+nguard)-((jd+jl)-(nyb+nguard+1))
                          jj  = jd + jl + 1
                       endif
                    endif
                 endif
              endif

              do j = ja,jb
                 ii = id
                 do i = ia,ib
                    if (k >= ks .and. k <= ks + kl) then
                       if (j >= jlo .and. j <= jhi) then
                          if (i >= is .and. i <= is + il) then

                             do ivar=1,ngcell_on_fc(2)
                                ivar_next = gcell_on_fc_pointer(2,ivar)

                                if (.not.l_f_to_c) then 
                                   if(j.eq.1+nguard0*k2d.and.iface.eq.4) then
                                      if (force_consistency) then
                                         recvyf(ivar_next,i,1,k) =  & 
                                              &              temprecv_buf(index0+ivar)
                                      endif
                                   elseif(j.eq.nyb+1+nguard0*k2d.and.iface.eq.3) then
                                      if (force_consistency) then
!!$                                         print*,'cp_rem:',recvyf(ivar_next,i,1+k2d,k),temprecv_buf(index0+ivar),ivar_next,i
                                         recvyf(ivar_next,i,1+k2d,k) =  & 
                                              &              temprecv_buf(index0+ivar)
                                      endif
                                   else
                                      facevary1(ivar_next,ii,jj,kk,idest) =  & 
                                           &              temprecv_buf(index0+ivar)
                                   endif
                                else
                                   facevary1_fl(ivar_next,ii,jj,kk) =  & 
                                        &              temprecv_buf(index0+ivar)
                                end if

                             enddo

                          end if
                       end if
                    end if
                    if (i >= is .and. i <= is + il) ii = ii + 1
                    index0 = index0+ngcell_on_fc(2)
                 enddo
                 !        if (j >= jlo .and. j <= jhi) jj = jj + 1
                 if (j >= jlo .and. j <= jhi) jj = jj + jstride
              enddo
              if (k >= ks .and. k <= ks + kl) kk = kk + 1
           enddo
        endif

       if(ndim.eq.3) then
        vtype = 4 
        call mpi_set_message_limits( & 
     &               dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &               ill,jll,kll)


        klo = ks1
        khi = ks1+kl+kp2
        if(iface.eq.5) khi = nzb+1+nguard0*k3d
        if(iface.eq.6) klo = 1+nguard0*npgs

        kk = kd1
        if(iface.eq.6) kk = kd1-1
        do k = ka,kb
        jj = jd
        jstride = 1

      js2 = js
      js1 = js+jl

      if(spherical_pm) then
      if(lsingular_line) then
      if(ipolar(1).eq.-1.and.jd.le.nguard) then
        jstride = -1
        if (no_permanent_guardcells) then
          js2 = jd + jl + (nguard - 2*(jd+jl)) +1
          js1 = jd      + (nguard - 2* jd    ) +1
          jj  = jd + jl
        else
          js2 = jd + jl + 2*(nguard - (jd+jl)) +1
          js1 = jd      + 2*(nguard -  jd    ) +1
          jj  = jd + jl
        endif
      elseif(ipolar(2).eq.+1.and.jd.gt.nyb+nguard) then
        jstride = -1
        if (no_permanent_guardcells) then
          js1 = nyb - ( jd    -(nyb+nguard+1))
          js2 = nyb - ((jd+jl)-(nyb+nguard+1))
          jj  = jd + jl
        else
          js1 = (nyb+nguard)-( jd    -(nyb+nguard+1))
          js2 = (nyb+nguard)-((jd+jl)-(nyb+nguard+1))
          jj  = jd + jl
        endif
      endif
      endif
      endif

        do j = ja,jb
        ii = id
        do i = ia,ib
           if (k >= klo .and. k <= khi) then
!           if (j >= js .and. j <= js + jl) then
           if (j >= js2 .and. j <= js1) then
           if (i >= is .and. i <= is + il) then

        do ivar=1,ngcell_on_fc(3)
          ivar_next = gcell_on_fc_pointer(3,ivar)

           if (.not.l_f_to_c) then
            if(k.eq.1+nguard0*k3d.and.iface.eq.6) then
               if (force_consistency) then
              recvzf(ivar_next,i,j,1) =  & 
     &              temprecv_buf(index0+ivar)
              endif
            elseif(k.eq.nzb+1+nguard0*k3d.and.iface.eq.5) then
               if (force_consistency) then
              recvzf(ivar_next,i,j,1+k3d) =  & 
     &              temprecv_buf(index0+ivar)
              endif
            else
              facevarz1(ivar_next,ii,jj,kk,idest) =  & 
     &              temprecv_buf(index0+ivar)
            endif
           else
            facevarz1_fl(ivar_next,ii,jj,kk) =  & 
     &              temprecv_buf(index0+ivar)
           end if

        enddo

           end if
           end if
           end if
           if (i >= is .and. i <= is + il) ii = ii + 1
          index0 = index0+ngcell_on_fc(3)
        enddo
!        if (j >= js .and. j <= js + jl) jj = jj + 1
        if (j >= js2 .and. j <= js1) jj = jj + jstride
        enddo
        if (k >= klo .and. k <= khi) kk = kk + 1
        enddo
       endif


       endif

     if (l_f_to_c) then
      if(idest.eq.1) then

        f2c_ind_facex(1,1,nblk_ind) = min( id1, & 
     &                                   f2c_ind_facex(1,1,nblk_ind))
        f2c_ind_facex(2,1,nblk_ind) = max( id1+il+ip2, & 
     &                                   f2c_ind_facex(2,1,nblk_ind))
        f2c_ind_facex(1,2,nblk_ind) = min( jd, & 
     &                                   f2c_ind_facex(1,2,nblk_ind))
        f2c_ind_facex(2,2,nblk_ind) = max( jd+jl, & 
     &                                   f2c_ind_facex(2,2,nblk_ind))
        f2c_ind_facex(1,3,nblk_ind) = min( kd, & 
     &                                   f2c_ind_facex(1,3,nblk_ind))
        f2c_ind_facex(2,3,nblk_ind) = max( kd+kl, & 
     &                                   f2c_ind_facex(2,3,nblk_ind))

        f2c_ind_facey(1,1,nblk_ind) = min( id, & 
     &                                   f2c_ind_facey(1,1,nblk_ind))
        f2c_ind_facey(2,1,nblk_ind) = max( id+il, & 
     &                                   f2c_ind_facey(2,1,nblk_ind))
        f2c_ind_facey(1,2,nblk_ind) = min( jd1, & 
     &                                   f2c_ind_facey(1,2,nblk_ind))
        f2c_ind_facey(2,2,nblk_ind) = max( jd1+jl+jp2, & 
     &                                   f2c_ind_facey(2,2,nblk_ind))
        f2c_ind_facey(1,3,nblk_ind) = min( kd, & 
     &                                   f2c_ind_facey(1,3,nblk_ind))
        f2c_ind_facey(2,3,nblk_ind) = max( kd+kl, & 
     &                                   f2c_ind_facey(2,3,nblk_ind))

        f2c_ind_facez(1,1,nblk_ind) = min( id, & 
     &                                   f2c_ind_facez(1,1,nblk_ind))
        f2c_ind_facez(2,1,nblk_ind) = max( id+il, & 
     &                                   f2c_ind_facez(2,1,nblk_ind))
        f2c_ind_facez(1,2,nblk_ind) = min( jd, & 
     &                                   f2c_ind_facez(1,2,nblk_ind))
        f2c_ind_facez(2,2,nblk_ind) = max( jd+jl, & 
     &                                   f2c_ind_facez(2,2,nblk_ind))
        f2c_ind_facez(1,3,nblk_ind) = min( kd1, & 
     &                                   f2c_ind_facez(1,3,nblk_ind))
        f2c_ind_facez(2,3,nblk_ind) = max( kd1+kl+kp2, & 
     &                                   f2c_ind_facez(2,3,nblk_ind))
      endif
      end if

! make sure srl shared faces end up with the mean value from each
! block^s data 
      if (force_consistency) then
      if(.not.l_f_to_c) then

       if(iface.eq.1) then
       do ivar=1,ngcell_on_fc(1)
         ivar_next = gcell_on_fc_pointer(1,ivar)

         facevarx1(ivar_next,1+nguard,1+nguard*k2d:nyb+nguard*k2d, & 
     &                       1+nguard*k3d:nzb+nguard*k3d,idest) = & 
     &   .5*( & 
     &   facevarx1(ivar_next,1+nguard,1+nguard*k2d:nyb+nguard*k2d, & 
     &                       1+nguard*k3d:nzb+nguard*k3d,idest) + & 
     &     recvxf(ivar_next,2,1+nguard0*k2d:nyb+nguard0*k2d, & 
     &                         1+nguard0*k3d:nzb+nguard0*k3d)  & 
     &      ) 
       enddo

       elseif(iface.eq.2) then
       do ivar=1,ngcell_on_fc(1)
         ivar_next = gcell_on_fc_pointer(1,ivar)

         facevarx1(ivar_next,nxb+1+nguard,1+nguard*k2d:nyb+nguard*k2d, & 
     &                       1+nguard*k3d:nzb+nguard*k3d,idest) = & 
     &   .5*( & 
     &     recvxf(ivar_next,1,1+nguard0*k2d:nyb+nguard0*k2d, & 
     &                       1+nguard0*k3d:nzb+nguard0*k3d)  & 
     &   + facevarx1(ivar_next,nxb+1+nguard,1+nguard*k2d:nyb+nguard*k2d, & 
     &                       1+nguard*k3d:nzb+nguard*k3d,idest)  & 
     &      )
       enddo

       elseif(iface.eq.3) then
        jbface = 2
      if(spherical_pm) then
        if(ipolar(1).eq.-1) jbface = 1
      endif
       do ivar=1,ngcell_on_fc(2)
         ivar_next = gcell_on_fc_pointer(2,ivar)

         facevary1(ivar_next,1+nguard:nxb+nguard,1+nguard*k2d, & 
     &                       1+nguard*k3d:nzb+nguard*k3d,idest) = & 
     &   .5*( & 
     &   facevary1(ivar_next,1+nguard:nxb+nguard,1+nguard*k2d, & 
     &                       1+nguard*k3d:nzb+nguard*k3d,idest) + & 
     &     recvyf(ivar_next,1+nguard0:nxb+nguard0,jbface, & 
     &                     1+nguard0*k3d:nzb+nguard0*k3d) )
       enddo

       elseif(iface.eq.4) then
        jbface = 1
      if(spherical_pm) then
        if(ipolar(2).eq.1) jbface = 2
      endif
       do ivar=1,ngcell_on_fc(2)
         ivar_next = gcell_on_fc_pointer(2,ivar)

         facevary1(ivar_next,1+nguard:nxb+nguard,nyb+k2d+nguard*k2d, & 
     &                       1+nguard*k3d:nzb+nguard*k3d,idest) = & 
     &   .5*( & 
     &     recvyf(ivar_next,1+nguard0:nxb+nguard0,jbface, & 
     &                     1+nguard0*k3d:nzb+nguard0*k3d)  & 
     &    + facevary1(ivar_next,1+nguard:nxb+nguard,nyb+k2d+nguard*k2d, & 
     &                       1+nguard*k3d:nzb+nguard*k3d,idest) & 
     &      )
       enddo

       elseif(iface.eq.5) then
       do ivar=1,ngcell_on_fc(3)
         ivar_next = gcell_on_fc_pointer(3,ivar)

         facevarz1(ivar_next,1+nguard:nxb+nguard, & 
     &             1+nguard*k2d:nyb+nguard*k2d,1+nguard*k3d,idest) = & 
     &   .5*( & 
     &   facevarz1(ivar_next,1+nguard:nxb+nguard, & 
     &             1+nguard*k2d:nyb+nguard*k2d,1+nguard*k3d,idest) + & 
     &     recvzf(ivar_next,1+nguard0:nxb+nguard0, & 
     &                     1+nguard0*k2d:nyb+nguard0*k2d,1+k3d) )
       enddo

       elseif(iface.eq.6) then
       do ivar=1,ngcell_on_fc(3)
         ivar_next = gcell_on_fc_pointer(3,ivar)

         facevarz1(ivar_next,1+nguard:nxb+nguard, & 
     &       1+nguard*k2d:nyb+nguard*k2d,nzb+k3d+nguard*k3d,idest) = & 
     &   .5*( & 
     &     recvzf(ivar_next,1+nguard0:nxb+nguard0, & 
     &                     1+nguard0*k2d:nyb+nguard0*k2d,1)  & 
     &   + facevarz1(ivar_next,1+nguard:nxb+nguard, & 
     &       1+nguard*k2d:nyb+nguard*k2d,nzb+k3d+nguard*k3d,idest) & 
     &       )

       enddo
       endif

      endif
      endif


      return
      end subroutine amr_1blk_fc_cp_remote
