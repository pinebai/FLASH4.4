      subroutine amr_edge_average_udt(mype)


! $RCSfile: amr_edge_average_udt.F90,v $
! $Revision: 1.1 $
! $Date: 2004/05/20 01:03:44 $


!------------------------------------------------------------------------
!
! This routine gets cell edge-based data at block boundaries from 
! neighbors who are parents of leaf blocks. 
!
! The data structure used to store and pass this data is defined
! in the include file 'block_boundary_data.fh' which can be included
! in 'physicaldata.fh'.
!
! This version is called when uniform timesteps are being used across
! the blocks in the computation.
!
!
! Written :     Peter MacNeice          August 1997
!------------------------------------------------------------------------

use physicaldata
      use tree
      implicit none



      integer mype

!------------------------------------
! local variables


      integer ng_off
      parameter(ng_off = nguard+iface_off)

      integer kup,klo,kup1
      parameter(klo  = 1+k3d*nguard)
      parameter(kup  = 1+k3d*(nzb+nguard-1))
      parameter(kup1 = k3d+nzb+k3d*nguard)

      integer remote_pe,remote_block
      integer cnodetype

      integer isg,jf

!------------------------------------


! all leaf blocks provide reduced boundary edge data to their parents
      call amr_restrict_edge_data(mype)

!      call shmem_barrier_all()


! cycle through the grid blocks on this processor
      if(lnblocks.gt.0) then
      do isg = 1,lnblocks

! Is this a leaf block and not at the original refinement level ?
      if(nodetype(isg).eq.1.and.lrefine(isg).gt.1) then


! Cycle over the blocks faces
        do jf = 1,nfaces

          remote_pe = neigh(2,jf,isg)
          remote_block  = neigh(1,jf,isg)

! Is the neighbor to this face a parent of a leaf block?
          cnodetype = 0
!          if(remote_block.gt.0) call shmem_integer_get(cnodetype,
!     .                         nodetype(remote_block),1,remote_pe)
          if(cnodetype.eq.2) then

! If yes then copy the appropriate layer from its boundary variable data 

            if(jf.eq.1) then

!              call shmem_real_get
!     . (recvarx1e(1,1,1,1),
!     .                       bedge_facex_y(1,1,1,1,remote_block),
!     .                       len_block_ex*nedges,remote_pe)
              bedge_facex_y(:,1,:,:,isg) = recvarx1e(:,2,:,:)

              if((ndim.eq.3).or.(l2p5d.eq.1)) then
!                call shmem_real_get
!     . (recvarx1e(1,1,1,1),
!     .                        bedge_facex_z(1,1,1,1,remote_block),
!     .                        len_block_ex*nedges,remote_pe)
                bedge_facex_z(:,1,:,:,isg) = recvarx1e(:,2,:,:)
              endif

! make common variables on an edge consistent
#if N_DIM >= 2
              bedge_facey_z(:,1+nguard,1,klo:kup,isg) = & 
     &               bedge_facex_z(:,1,1+nguard,klo:kup,isg)
       
              bedge_facey_z(:,1+nguard,2,klo:kup,isg) = & 
     &               bedge_facex_z(:,1,1+nyb+nguard,klo:kup,isg)

#endif
#if N_DIM == 3
                bedge_facez_y(:,1+nguard,1+nguard:nyb+nguard,1,isg) & 
     &             = bedge_facex_y(:,1,1+nguard:nyb+nguard,klo,isg)

                bedge_facez_y(:,1+nguard,1+nguard:nyb+nguard,2,isg) & 
     &             = bedge_facex_y(:,1,1+nguard:nyb+nguard,kup1,isg)
#endif


            elseif(jf.eq.2) then

!              call shmem_real_get
!     . (recvarx1e(1,1,1,1),
!     .               bedge_facex_y(1,1,1,1,remote_block),
!     .               len_block_ex*nedges,remote_pe)
              bedge_facex_y(:,2,:,:,isg)=recvarx1e(:,1,:,:)

              if((ndim.eq.3).or.(l2p5d.eq.1)) then
!                call shmem_real_get
!     . (recvarx1e(1,1,1,1),
!     .               bedge_facex_z(1,1,1,1,remote_block),
!     .               len_block_ex*nedges,remote_pe)
                bedge_facex_z(:,2,:,:,isg)=recvarx1e(:,1,:,:)
              endif

! make common variables on an edge consistent
#if N_DIM >= 2
              bedge_facey_z(:,1+nxb+nguard,1,klo:kup,isg) = & 
     &            bedge_facex_z(:,2,1+nguard,klo:kup,isg)
                
              bedge_facey_z(:,1+nxb+nguard,2,klo:kup,isg) = & 
     &            bedge_facex_z(:,2,1+nyb+nguard,klo:kup,isg)

#endif
#if N_DIM == 3
                bedge_facez_y(:,1+nxb+nguard,1+nguard:nyb+nguard,1,isg)= & 
     &            bedge_facex_y(:,2,1+nguard:nyb+nguard,klo,isg)

                bedge_facez_y(:,1+nxb+nguard,1+nguard:nyb+nguard,2,isg)= & 
     &            bedge_facex_y(:,2,1+nguard:nyb+nguard,kup1,isg)
#endif

            elseif(jf.eq.3) then

              if((ndim.eq.3).or.(l2p5d.eq.1)) then
!                call shmem_real_get
!     . (recvary1e(1,1,1,1),
!     .               bedge_facey_z(1,1,1,1,remote_block),
!     .               len_block_ey*nedges,remote_pe)
                bedge_facey_z(:,:,1,:,isg) = recvary1e(:,:,2,:)
              endif

!              call shmem_real_get
!     . (recvary1e(1,1,1,1),
!     .               bedge_facey_x(1,1,1,1,remote_block),
!     .               len_block_ey*nedges,remote_pe)
              bedge_facey_x(:,:,1,:,isg) = recvary1e(:,:,2,:)


! make common variables on an edge consistent
#if N_DIM >= 2
              bedge_facex_z(:,1,1+nguard,klo:kup,isg) = & 
     &          bedge_facey_z(:,1+nguard,1,klo:kup,isg)

              bedge_facex_z(:,2,1+nguard,klo:kup,isg) = & 
     &          bedge_facey_z(:,1+nxb+nguard,1,klo:kup,isg)

#endif
#if N_DIM == 3
                bedge_facez_x(:,1+nguard:nxb+nguard,1+nguard,1,isg)= & 
     &          bedge_facey_x(:,1+nguard:nxb+nguard,1,klo,isg)

                bedge_facez_x(:,1+nguard:nxb+nguard,1+nguard,2,isg)= & 
     &          bedge_facey_x(:,1+nguard:nxb+nguard,1,kup1,isg)
#endif

              elseif(jf.eq.4) then

                if((ndim.eq.3).or.(l2p5d.eq.1)) then
!                  call shmem_real_get
!     . (recvary1e(1,1,1,1),
!     .               bedge_facey_z(1,1,1,1,remote_block),
!     .               len_block_ey*nedges,remote_pe)
                  bedge_facey_z(:,:,2,:,isg) = recvary1e(:,:,1,:)
                endif

!                call shmem_real_get
!     . (recvary1e(1,1,1,1),
!     .               bedge_facey_x(1,1,1,1,remote_block),
!     .               len_block_ey*nedges,remote_pe)
                bedge_facey_x(:,:,2,:,isg) = recvary1e(:,:,1,:)


! make common variables on an edge consistent
#if N_DIM >= 2
                bedge_facex_z(:,1,1+nyb+nguard,klo:kup,isg) = & 
     &              bedge_facey_z(:,1+nguard,2,klo:kup,isg)

                bedge_facex_z(:,2,1+nyb+nguard,klo:kup,isg) = & 
     &              bedge_facey_z(:,1+nxb+nguard,2,klo:kup,isg)

#endif
#if N_DIM == 3
                bedge_facez_x(:,1+nguard:nxb+nguard,1+nyb+nguard,1,isg)= & 
     &              bedge_facey_x(:,1+nguard:nxb+nguard,2,klo,isg)

                bedge_facez_x(:,1+nguard:nxb+nguard,1+nyb+nguard,2,isg)= & 
     &              bedge_facey_x(:,1+nguard:nxb+nguard,2,kup1,isg)
#endif


              elseif(jf.eq.5) then

!                call shmem_real_get
!     . (recvarz1e(1,1,1,1),
!     .               bedge_facez_x(1,1,1,1,remote_block),
!     .               len_block_ez*nedges,remote_pe)
                bedge_facez_x(:,:,:,1,isg) = recvarz1e(:,:,:,2)
!                call shmem_real_get
!     . (recvarz1e(1,1,1,1),
!     .               bedge_facez_y(1,1,1,1,remote_block),
!     .               len_block_ez*nedges,remote_pe)
                bedge_facez_y(:,:,:,1,isg) = recvarz1e(:,:,:,2) 

! make common variables on an edge consistent
#if N_DIM >= 2
                bedge_facey_x(:,1+nguard:nxb+nguard,1,klo,isg)= & 
     &          bedge_facez_x(:,1+nguard:nxb+nguard,1+nguard,1,isg)

                bedge_facey_x(:,1+nguard:nxb+nguard,2,klo,isg)= & 
     &          bedge_facez_x(:,1+nguard:nxb+nguard,1+nyb+nguard,1,isg)

                bedge_facex_y(:,1,1+nguard:nyb+nguard,klo,isg)= & 
     &          bedge_facez_y(:,1+nguard,1+nguard:nyb+nguard,1,isg)

                bedge_facex_y(:,2,1+nguard:nyb+nguard,klo,isg)= & 
     &          bedge_facez_y(:,1+nxb+nguard,1+nguard:nyb+nguard,1,isg)
#endif

              elseif(jf.eq.6) then 

!                call shmem_real_get
!     . (recvarz1e(1,1,1,1),
!     .               bedge_facez_x(1,1,1,1,remote_block),
!     .               len_block_ez*nedges,remote_pe)
                bedge_facez_x(:,:,:,2,isg) = recvarz1e(:,:,:,1)
!                call shmem_real_get
!     . (recvarz1e(1,1,1,1),
!     .               bedge_facez_y(1,1,1,1,remote_block),
!     .               len_block_ez*nedges,remote_pe)
                bedge_facez_y(:,:,:,2,isg) = recvarz1e(:,:,:,1)

! make common variables on an edge consistent
#if N_DIM >= 2
                bedge_facey_x(:,1+nguard:nxb+nguard,1,kup1,isg)= & 
     &          bedge_facez_x(:,1+nguard:nxb+nguard,1+nguard,2,isg)

                bedge_facey_x(:,1+nguard:nxb+nguard,2,kup1,isg)= & 
     &          bedge_facez_x(:,1+nguard:nxb+nguard,1+nyb+nguard,2,isg)

                bedge_facex_y(:,1,1+nguard:nyb+nguard,kup1,isg)= & 
     &          bedge_facez_y(:,1+nguard,1+nguard:nyb+nguard,2,isg)

                bedge_facex_y(:,2,1+nguard:nyb+nguard,kup1,isg)= & 
     &          bedge_facez_y(:,1+nxb+nguard,1+nguard:nyb+nguard,2,isg)
#endif


              endif


          endif

        enddo

      endif
      enddo
      endif

!      call shmem_barrier_all()

      return
      end
