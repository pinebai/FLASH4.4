      subroutine amr_prolong_fc_divbconsist(mype)


! $RCSfile: amr_prolong_fc_divbconsist.F90,v $
! $Revision: 1.1 $
! $Date: 2004/05/20 01:03:44 $


!------------------------------------------------------------------------
!
! Like the routine amr_prolong_fc_consist this routine checks for
! existing neighbor to newly created child blocks, and where found,
! uses facevar data from the existing neighbor at the common block
! boundary, in place of interpolation from the new childs parent.
! In addition this routine makesthis change while adjusting values
! immediately interior to the new face to guarantee that div B is
! kept at zero.
!
! Written :     Peter MacNeice          April 1998
!------------------------------------------------------------------------

use physicaldata
      use tree
       implicit none


!------------------------------------
! local variables

      integer mype
      integer isw,ipe,ibl,ii,jj,kk,idv,i,j,k
      integer isg,jf
      integer idest
      integer i_dest,i_source,j_dest,j_source,k_dest,k_source

      real bsum
      real divbmax,divb

      integer remote_pe,remote_block
      logical cnewchild
      save    cnewchild

! local arrays
      real recvx(nbndvar,il_bnd:iu_bnd+1,jl_bnd:ju_bnd, & 
     &       kl_bnd:ku_bnd)
      real recvy(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d, & 
     &       kl_bnd:ku_bnd)
      real recvz(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd, & 
     &       kl_bnd:ku_bnd+k3d)
      save recvx,recvy,recvz

      real tempx(jl_bnd:ju_bnd,kl_bnd:ku_bnd)
      real tempy(il_bnd:iu_bnd,kl_bnd:ku_bnd)
      real tempz(il_bnd:iu_bnd,jl_bnd:ju_bnd)

      real e1,e2,ea,eb
      real b11,b12,b21,b22

!------------------------------------

!      call shmem_barrier_all()

      isw = 0

      ipe = 152
      ibl = 1
      ii = 3
      jj = 10
      kk = 10

!
! select component of facevar to which div correction is to be applied.
      idv = 1

! cycle through the grid blocks on this processor
      if(lnblocks.gt.0) then
      do isg = 1,lnblocks

! Is this a newly created leaf block ?
      if(nodetype(isg).eq.1.and.newchild(isg)) then


! Cycle over the blocks faces
       do jf = 1,nfaces

          remote_pe = neigh(2,jf,isg)
          remote_block  = neigh(1,jf,isg)

! Is the neighbor to this face a pre-existing block?
          cnewchild = .true.
!          if(remote_block.gt.0) call shmem_logical_get(cnewchild,
!     .                      newchild(remote_block),1,remote_pe)

          if(.not.cnewchild) then

! If the neighbor block is pre-existing then get its facevar data on the
! shared block boundary.

            idest = isg

            if(jf.eq.1) then
               i_dest   = nguard + 1 + iface_off
               i_source = nxb+nguard + 1 -gc_off_x + iface_off
               tempx(:,:) = facevarx(idv,i_dest,:,:,idest)

!           if(mype.eq.ipe.and.isg.eq.ibl) then
!           write(*,*) 'diag: old bx '
!     .     ,facevarx(1,3,9,9,idest),facevarx(1,3,10,9,idest)
!     .     ,facevarx(1,3,9,10,idest),facevarx(1,3,10,10,idest)
!           endif
!               call shmem_real_get (
!     .                    recvx,facevarx(1,1,1,1,remote_block),
!     .                    len_blockfx*nbndvar,remote_pe)
               facevarx(idv,i_dest,:,:,idest)= recvx(idv,i_source,:,:)
!           if(mype.eq.ipe.and.isg.eq.ibl) then
!           write(*,*) 'diag: new bx '
!     .     ,facevarx(1,3,9,9,idest),facevarx(1,3,10,9,idest)
!     .     ,facevarx(1,3,9,10,idest),facevarx(1,3,10,10,idest)
!           endif

               do k=nguard+1,nguard + nzb-1,2
               do j=nguard+1,nguard + nyb-1,2
                b11 = facevarx(idv,nguard+1,j,k,idest) & 
     &                                      - tempx(j,k)
                b21 = facevarx(idv,nguard+1,j,k+1,idest) & 
     &                                      - tempx(j,k+1)
                b12 = facevarx(idv,nguard+1,j+1,k,idest) & 
     &                                      - tempx(j+1,k)
                b22 = facevarx(idv,nguard+1,j+1,k+1,idest) & 
     &                                      - tempx(j+1,k+1)
                bsum = b11+b12+b21+b22
!           if(mype.eq.ipe.and.isg.eq.ibl.and.j.eq.9.and.k.eq.9) then
!           write(*,*) 'diag: db ',b11,b12,b21,b22,bsum
!           isw = 1
!           endif

!                if(bsum.gt.1.e-5) write(*,*) 'bsum : ',
!     .                                    mype,isg,jf,j,k,bsum
                call compute_evalues(b11,b12,b21,b22,ea,eb,e1,e2,isw)

!           if(mype.eq.ipe.and.isg.eq.ibl.and.j.eq.9.and.k.eq.9) then
!           write(*,*) 'diag: old by/z '
!     .     ,facevary(1,3,10,9,idest),facevary(1,3,10,10,idest)
!     .     ,facevarz(1,3,9,10,idest),facevarz(1,3,10,10,idest)
!          endif
                 facevary(idv,nguard+1,j+1,k,idest) = & 
     &             facevary(idv,nguard+1,j+1,k,idest)-e1
                 facevary(idv,nguard+1,j+1,k+1,idest) = & 
     &             facevary(idv,nguard+1,j+1,k+1,idest)-e2
                 facevarz(idv,nguard+1,j,k+1,idest) = & 
     &             facevarz(idv,nguard+1,j,k+1,idest)+ea
                 facevarz(idv,nguard+1,j+1,k+1,idest) = & 
     &             facevarz(idv,nguard+1,j+1,k+1,idest)+eb
!           if(mype.eq.ipe.and.isg.eq.ibl.and.j.eq.3.and.k.eq.5) then
!           write(*,*) 'diag: es ',ea,eb,e1,e2
!           write(*,*) 'diag: modif by/z '
!     .     ,facevary(1,3,10,9,idest),facevary(1,3,10,10,idest)
!     .     ,facevarz(1,3,9,10,idest),facevarz(1,3,10,10,idest)
!           endif

               enddo
               enddo

!           if(mype.eq.ipe.and.isg.eq.ibl) then
!           write(*,*) 'diag: other bx '
!     .     ,facevarx(1,4,9,9,idest),facevarx(1,4,10,9,idest)
!     .     ,facevarx(1,4,9,10,idest),facevarx(1,4,10,10,idest)
!           write(*,*) 'diag: other by low side '
!     .     ,facevary(1,3,9,9,idest),facevary(1,3,9,10,idest)
!           write(*,*) 'diag: other by hi side '
!     .     ,facevary(1,3,11,9,idest),facevary(1,3,11,10,idest)
!           write(*,*) 'diag: other bz low side '
!     .     ,facevarz(1,3,9,9,idest),facevarz(1,3,10,9,idest)
!           write(*,*) 'diag: other bz hi side '
!     .     ,facevarz(1,3,9,11,idest),facevarz(1,3,10,11,idest)
!             i=3
!             do k=9,10
!             do j=9,10
!                 divb = (
!     .            facevarx(1,i+1,j,k,isg) - facevarx(1,i,j,k,isg)
!     .            + facevary(1,i,j+1,k,isg) - facevary(1,i,j,k,isg)
!     .            + facevarz(1,i,j,k+1,isg) - facevarz(1,i,j,k,isg))
!             write(*,*) 'diag : divB ',i,j,k,divb
!             enddo
!             enddo
!           endif

            elseif(jf.eq.2) then

               i_dest   = nxb+1+nguard + iface_off
               i_source = 1+nguard+gc_off_x + iface_off
               tempx(:,:) = facevarx(idv,i_dest,:,:,idest)
!               call shmem_real_get (
!     .                    recvx,facevarx(1,1,1,1,remote_block),
!     .                    len_blockfx*nbndvar,remote_pe)
               facevarx(idv,i_dest,:,:,idest)= recvx(idv,i_source,:,:)

               do k=nguard+1,nguard + nzb-1,2
               do j=nguard+1,nguard + nyb-1,2
                b11 = facevarx(idv,nguard+nxb+1,j,k,idest) & 
     &                                              - tempx(j,k)
                b21 = facevarx(idv,nguard+nxb+1,j,k+1,idest) & 
     &                                              - tempx(j,k+1)
                b12 = facevarx(idv,nguard+nxb+1,j+1,k,idest) & 
     &                                              - tempx(j+1,k)
                b22 = facevarx(idv,nguard+nxb+1,j+1,k+1,idest) & 
     &                                              - tempx(j+1,k+1)
                bsum = b11+b12+b21+b22
!                if(bsum.gt.1.e-5) write(*,*) 'bsum : ',
!     .                                    mype,isg,jf,j,k,bsum
                call compute_evalues(b11,b12,b21,b22,ea,eb,e1,e2,isw)

                facevary(idv,nguard+nxb,j+1,k,idest) = & 
     &             facevary(idv,nguard+nxb,j+1,k,idest)+e1
                facevary(idv,nguard+nxb,j+1,k+1,idest) = & 
     &             facevary(idv,nguard+nxb,j+1,k+1,idest)+e2
                facevarz(idv,nguard+nxb,j,k+1,idest) = & 
     &             facevarz(idv,nguard+nxb,j,k+1,idest)-ea
                facevarz(idv,nguard+nxb,j+1,k+1,idest) = & 
     &             facevarz(idv,nguard+nxb,j+1,k+1,idest)-eb
               enddo
               enddo


#if N_DIM >= 2
            elseif(jf.eq.3) then

               j_dest   = nguard + 1 + iface_off
               j_source = nyb+nguard + 1 -gc_off_y + iface_off
               tempy(:,:) = facevary(idv,:,j_dest,:,idest)
!               call shmem_real_get (
!     .                    recvy,facevary(1,1,1,1,remote_block),
!     .                    len_blockfy*nbndvar,remote_pe)
               facevary(idv,:,j_dest,:,idest)=recvy(idv,:,j_source,:)

               do k=nguard+1,nguard + nzb-1,2
               do i=nguard+1,nguard + nxb-1,2
                b11 = facevary(idv,i,nguard+1,k,idest) - tempy(i,k)
                b21 = facevary(idv,i,nguard+1,k+1,idest) & 
     &                                             - tempy(i,k+1)
                b12 = facevary(idv,i+1,nguard+1,k,idest) & 
     &                                             - tempy(i+1,k)
                b22 = facevary(idv,i+1,nguard+1,k+1,idest) & 
     &                                             - tempy(i+1,k+1)
                bsum = b11+b12+b21+b22
!                if(bsum.gt.1.e-5) write(*,*) 'bsum : ',
!     .                                    mype,isg,jf,i,k,bsum
                call compute_evalues(b11,b12,b21,b22,ea,eb,e1,e2,isw)

                facevarx(idv,i+1,nguard+1,k,idest) = & 
     &             facevarx(idv,i+1,nguard+1,k,idest)-e1
                facevarx(idv,i+1,nguard+1,k+1,idest) = & 
     &             facevarx(idv,i+1,nguard+1,k+1,idest)-e2
                facevarz(idv,i,nguard+1,k+1,idest) = & 
     &             facevarz(idv,i,nguard+1,k+1,idest)+ea
                facevarz(idv,i+1,nguard+1,k+1,idest) = & 
     &             facevarz(idv,i+1,nguard+1,k+1,idest)+eb
               enddo
               enddo

            elseif(jf.eq.4) then

               j_dest   = nyb+1+nguard + iface_off
               j_source = 1+nguard+gc_off_y + iface_off
               tempy(:,:) = facevary(idv,:,j_dest,:,idest)
!               call shmem_real_get (
!     .                    recvy,facevary(1,1,1,1,remote_block),
!     .                    len_blockfy*nbndvar,remote_pe)
               facevary(idv,:,j_dest,:,idest)=recvy(idv,:,j_source,:)

               do k=nguard+1,nguard + nzb-1,2
               do i=nguard+1,nguard + nxb-1,2
                b11 = facevary(idv,i,nguard+nyb+1,k,idest) & 
     &                                        - tempy(i,k)
                b21 = facevary(idv,i,nguard+nyb+1,k+1,idest) & 
     &                                        - tempy(i,k+1)
                b12 = facevary(idv,i+1,nguard+nyb+1,k,idest) & 
     &                                        - tempy(i+1,k)
                b22 = facevary(idv,i+1,nguard+nyb+1,k+1,idest) & 
     &                                        - tempy(i+1,k+1)
                bsum = b11+b12+b21+b22
!                if(bsum.gt.1.e-5) write(*,*) 'bsum : ',
!     .                                    mype,isg,jf,i,k,bsum
                call compute_evalues(b11,b12,b21,b22,ea,eb,e1,e2,isw)

                facevarx(idv,i+1,nguard+nyb,k,idest) = & 
     &             facevarx(idv,i+1,nguard+nyb,k,idest)+e1
                facevarx(idv,i+1,nguard+nyb,k+1,idest) = & 
     &             facevarx(idv,i+1,nguard+nyb,k+1,idest)+e2
                facevarz(idv,i,nguard+nyb,k+1,idest) = & 
     &             facevarz(idv,i,nguard+nyb,k+1,idest)-ea
                facevarz(idv,i+1,nguard+nyb,k+1,idest) = & 
     &             facevarz(idv,i+1,nguard+nyb,k+1,idest)-eb
               enddo
               enddo

#endif
#if N_DIM == 3
            elseif(jf.eq.5) then

               k_dest   = nguard + 1 + iface_off
               k_source = nzb+nguard + 1 -gc_off_z + iface_off
               tempz(:,:) = facevarz(idv,:,:,k_dest,idest)
!               call shmem_real_get (
!     .                    recvz,facevarz(1,1,1,1,remote_block),
!     .                    len_blockfz*nbndvar,remote_pe)
               facevarz(idv,:,:,k_dest,idest)=recvz(idv,:,:,k_source)

               do j=nguard+1,nguard + nyb-1,2
               do i=nguard+1,nguard + nxb-1,2
                b11 = facevarz(idv,i,j,nguard+1,idest) - tempz(i,j)
                b21 = facevarz(idv,i,j+1,nguard+1,idest) & 
     &                                          - tempz(i,j+1)
                b12 = facevarz(idv,i+1,j,nguard+1,idest) & 
     &                                          - tempz(i+1,j)
                b22 = facevarz(idv,i+1,j+1,nguard+1,idest) & 
     &                                          - tempz(i+1,j+1)
                bsum = b11+b12+b21+b22
!                if(bsum.gt.1.e-5) write(*,*) 'bsum : ',
!     .                                    mype,isg,jf,i,j,bsum
                call compute_evalues(b11,b12,b21,b22,ea,eb,e1,e2,isw)

                facevarx(idv,i+1,j,nguard+1,idest) = & 
     &             facevarx(idv,i+1,j,nguard+1,idest)-e1
                facevarx(idv,i+1,j+1,nguard+1,idest) = & 
     &             facevarx(idv,i+1,j+1,nguard+1,idest)-e2
                facevary(idv,i,j+1,nguard+1,idest) = & 
     &             facevary(idv,i,j+1,nguard+1,idest)+ea
                facevary(idv,i+1,j+1,nguard+1,idest) = & 
     &             facevary(idv,i+1,j+1,nguard+1,idest)+eb
               enddo
               enddo

            elseif(jf.eq.6) then

               k_dest   = nzb+1+nguard + iface_off
               k_source = 1+nguard+gc_off_z + iface_off
               tempz(:,:) = facevarz(idv,:,:,k_dest,idest)
!               call shmem_real_get (
!     .                    recvz,facevarz(1,1,1,1,remote_block),
!     .                    len_blockfz*nbndvar,remote_pe)
               facevarz(idv,:,:,k_dest,idest)=recvz(idv,:,:,k_source)

               do j=nguard+1,nguard + nyb -1,2
               do i=nguard+1,nguard + nxb -1,2
                b11 = facevarz(idv,i,j,nguard+nzb+1,idest) & 
     &                                            - tempz(i,j)
                b21 = facevarz(idv,i,j+1,nguard+nzb+1,idest) & 
     &                                            - tempz(i,j+1)
                b12 = facevarz(idv,i+1,j,nguard+nzb+1,idest) & 
     &                                            - tempz(i+1,j)
                b22 = facevarz(idv,i+1,j+1,nguard+nzb+1,idest) & 
     &                                            - tempz(i+1,j+1)
                bsum = b11+b12+b21+b22
!                if(bsum.gt.1.e-5) write(*,*) 'bsum : ',
!     .                                    mype,isg,jf,i,j,bsum
                call compute_evalues(b11,b12,b21,b22,ea,eb,e1,e2,isw)

                facevarx(idv,i+1,j,nguard+nzb,idest) = & 
     &             facevarx(idv,i+1,j,nguard+nzb,idest)+e1
                facevarx(idv,i+1,j+1,nguard+nzb,idest) = & 
     &             facevarx(idv,i+1,j+1,nguard+nzb,idest)+e2
                facevary(idv,i,j+1,nguard+nzb,idest) = & 
     &             facevary(idv,i,j+1,nguard+nzb,idest)-ea
                facevary(idv,i+1,j+1,nguard+nzb,idest) = & 
     &             facevary(idv,i+1,j+1,nguard+nzb,idest)-eb
               enddo
               enddo


#endif
             endif

          endif

       enddo

      endif
      enddo
      endif

!      call shmem_barrier_all()


!
! div B check

      if(lnblocks.gt.0) then
      do isg = 1,lnblocks

! Is this a newly created leaf block ?
      if(nodetype(isg).eq.1) then

               divbmax = 0.
               do k=nguard+1,nguard + nzb
               do j=nguard+1,nguard + nyb
               do i=nguard+1,nguard + nxb
                 divb = ( & 
     &            facevarx(1,i+1,j,k,isg) - facevarx(1,i,j,k,isg) & 
     &            + facevary(1,i,j+1,k,isg) - facevary(1,i,j,k,isg) & 
     &            + facevarz(1,i,j,k+1,isg) - facevarz(1,i,j,k,isg))
                divbmax = max(divbmax,abs(divb))
               enddo
               enddo
               enddo

!               write(*,*) 'Max div B on block ',isg,' is ',divbmax
      endif
      enddo
      endif

!      call shmem_barrier_all()

      return
      end



      subroutine  compute_evalues(b11,b12,b21,b22,ea,eb,e1,e2,isw)

!
! Computes virtual electric field values required to produce
! changes in facevar's in a manner which will preserve div = 0
! constraints. The b11, b12 etc input arguments specify the
! required change in the vector component on the chosen face.
! The ea,eb,e1,e2 arguments return the virtual electric field
! values required to achieve this adjustment.
!
!
! The relationship between b's and e's is
!
!               -------------------------------------
!              |                  |                  |
!              |                  |                  |
!              |                  ^                  |
!              |       b21        e2     b22         |
!              |                  |                  |
!              |                  |                  |
!              |------ ea ->------|----- eb ->-------|
!              |                  |                  |
!              |                  ^                  |
!              |       b11        e1     b12         |
!              |                  |                  |
!              |                  |                  |
!              |                  |                  |
!               -------------------------------------
!
! These electric fields must be applied to adjust the appropriate
! components of the vector field in the planes perpendicular to
! the chosen face, immediately inside this face.


use physicaldata


! local variables

      real e1,e2,ea,eb
      real b11,b12,b21,b22
      integer isw
      

      ea = ( 2.*b11 +    b12 -    b21          )*.25
      eb = (    b11 + 2.*b12             - b22 )*.25
      e1 = (-2.*b11 +    b12 -    b21          )*.25
      e2 = (   -b11          - 2.*b21 + b22    )*.25

      if(isw.eq.1) then
       write(*,*) 'compute : b ',b11,b12,b21,b22
       write(*,*) 'compute : e ',ea,eb,e1,e2
       isw = 0
      endif

      return
      end
