      subroutine amr_edge_diagonal_check(mype)


! $RCSfile: amr_edge_diagonal_check.F90,v $
! $Revision: 1.1 $
! $Date: 2004/05/20 01:03:44 $


!------------------------------------------------------------------------
!
! This routine checks to see if the diagonal block between two
! leaf-neighbors at the same refinement level as the current block,
! is refined. If it is then the edge-based variables along the edge
! shared with that diagonal block is given the edge values
! form the refined diagonal block, to insure conservation properties.
!
!
! Written :     Peter MacNeice          October 1997
!------------------------------------------------------------------------

use physicaldata
      use tree
      implicit none


      integer mype

!------------------------------------
! local variables

      integer klo,kup
      parameter(klo=1+nguard*k3d,kup=klo*k3d+nzb)
      integer jlo,jup
      parameter(jlo=1+nguard,jup=jlo+nyb)
      integer ilo,iup
      parameter(ilo=1+nguard,iup=ilo+nxb)

      real receive(nedges,maxdim+2*nguard)

      integer remote_pe,remote_block
      integer neigh_ref_lev(mfaces,maxblocks),tneigh_ref_lev(mfaces)
      integer leaf_level(maxblocks)
      integer mark_edge(12,maxblocks),mark_face(6)

      integer lb,jf,l_ref,ie,i,j,k

      save neigh_ref_lev,tneigh_ref_lev,leaf_level
      save receive

!------------------------------------


! Initialize arrays storing refinement levels of neighboring leaf blocks,
! and marking edges for diagonal patching.
      neigh_ref_lev(:,:) = 0
      mark_edge(:,:) = 0


! Loop over the blocks on this processor. Record the leaf refinement levels
! associated with leaf and leaf-parent blocks.
      if(lnblocks.gt.0) then
      do lb=1,lnblocks
       leaf_level(lb) = lrefine(lb)
       if(nodetype(lb).eq.2) leaf_level(lb) = lrefine(lb)+1
      enddo
      endif


!        call shmem_barrier_all()

! Loop over the blocks on this processor.
      if(lnblocks.gt.0) then
      do lb=1,lnblocks

! Is this a leaf block ?
      if(nodetype(lb).eq.1) then


! Cycle over the blocks faces
       do jf = 1,nfaces

                remote_pe    = neigh(2,jf,lb)
                remote_block = neigh(1,jf,lb)
! Find the neighboring leaf refinement level.
!       if(remote_block.gt.0) call shmem_integer_get(
!     .       neigh_ref_lev(jf,lb),
!     .       leaf_level(remote_block),1,remote_pe)
       enddo

      l_ref = lrefine(lb)

! Test refinement levels of the leaf blocks facing the current block. If
! an edge is bounded by faces with the same refinement level as the current
! block, mark that edge as a possible candidate for a diagonal patch.

      if( ( neigh_ref_lev(1,lb).eq.l_ref) .and. & 
     &       (neigh_ref_lev(3,lb).eq.l_ref)) mark_edge(1,lb)=1
      if( ( neigh_ref_lev(1,lb).eq.l_ref) .and. & 
     &       (neigh_ref_lev(4,lb).eq.l_ref)) mark_edge(2,lb)=1
      if( ( neigh_ref_lev(2,lb).eq.l_ref) .and. & 
     &       (neigh_ref_lev(3,lb).eq.l_ref)) mark_edge(3,lb)=1
      if( ( neigh_ref_lev(2,lb).eq.l_ref) .and. & 
     &       (neigh_ref_lev(4,lb).eq.l_ref)) mark_edge(4,lb)=1

      if(ndim.eq.3) then

      if( ( neigh_ref_lev(3,lb).eq.l_ref) .and. & 
     &       (neigh_ref_lev(5,lb).eq.l_ref)) mark_edge(5,lb)=1
      if( ( neigh_ref_lev(4,lb).eq.l_ref) .and. & 
     &       (neigh_ref_lev(5,lb).eq.l_ref)) mark_edge(6,lb)=1
      if( ( neigh_ref_lev(3,lb).eq.l_ref) .and. & 
     &       (neigh_ref_lev(6,lb).eq.l_ref)) mark_edge(7,lb)=1
      if( ( neigh_ref_lev(4,lb).eq.l_ref) .and. & 
     &       (neigh_ref_lev(6,lb).eq.l_ref)) mark_edge(8,lb)=1
      if( ( neigh_ref_lev(1,lb).eq.l_ref) .and. & 
     &       (neigh_ref_lev(5,lb).eq.l_ref)) mark_edge(9,lb)=1
      if( ( neigh_ref_lev(2,lb).eq.l_ref) .and. & 
     &       (neigh_ref_lev(5,lb).eq.l_ref)) mark_edge(11,lb)=1
      if( ( neigh_ref_lev(1,lb).eq.l_ref) .and. & 
     &       (neigh_ref_lev(6,lb).eq.l_ref)) mark_edge(10,lb)=1
      if( ( neigh_ref_lev(2,lb).eq.l_ref) .and. & 
     &       (neigh_ref_lev(6,lb).eq.l_ref)) mark_edge(12,lb)=1

      endif


      endif

      enddo
      endif


!        call shmem_barrier_all()

! Loop over the blocks on this processor.
      if(lnblocks.gt.0) then
      do lb=1,lnblocks

       mark_face(:) = 0

! Is this a leaf block ?
      if(nodetype(lb).eq.1) then

! Set a marker on any face which has an edge already marked.
       mark_face(1) = mark_edge(1,lb)+mark_edge(2,lb) & 
     &       +mark_edge(9,lb)+mark_edge(10,lb)
       mark_face(2) = mark_edge(3,lb)+mark_edge(4,lb) & 
     &       +mark_edge(11,lb)+mark_edge(12,lb)
       mark_face(3) = mark_edge(1,lb)+mark_edge(3,lb) & 
     &       +mark_edge(5,lb)+mark_edge(7,lb)
       mark_face(4) = mark_edge(2,lb)+mark_edge(4,lb) & 
     &       +mark_edge(6,lb)+mark_edge(8,lb)
       if(ndim.eq.3) then
       mark_face(5) = mark_edge(5,lb)+mark_edge(6,lb) & 
     &       +mark_edge(9,lb)+mark_edge(11,lb)
       mark_face(6) = mark_edge(7,lb)+mark_edge(8,lb) & 
     &       +mark_edge(10,lb)+mark_edge(12,lb)
       endif

! Loop over all the faces of the current block.
       do jf = 1,nfaces

! Does this face have any edges marked?
       if(mark_face(jf).gt.0) then

       l_ref = lrefine(lb)+1

! Get the neigh_ref_lev array from this neighbor.
                remote_pe    = neigh(2,jf,lb)
                remote_block = neigh(1,jf,lb)
       tneigh_ref_lev(:)=0
!       if(remote_block.gt.0) call shmem_integer_get( 
!     .       tneigh_ref_lev,neigh_ref_lev(1,remote_block),
!     .       mfaces,remote_pe)

       if(jf.eq.1) then
       if(tneigh_ref_lev(3).ne.l_ref) mark_edge(1,lb)=0
       if(tneigh_ref_lev(4).ne.l_ref) mark_edge(2,lb)=0
       if(ndim.eq.3) then
       if(tneigh_ref_lev(5).ne.l_ref) mark_edge(9,lb)=0
       if(tneigh_ref_lev(6).ne.l_ref) mark_edge(10,lb)=0
       endif
       elseif(jf.eq.2) then
       if(tneigh_ref_lev(3).ne.l_ref) mark_edge(3,lb)=0
       if(tneigh_ref_lev(4).ne.l_ref) mark_edge(4,lb)=0
       if(ndim.eq.3) then
       if(tneigh_ref_lev(5).ne.l_ref) mark_edge(11,lb)=0
       if(tneigh_ref_lev(6).ne.l_ref) mark_edge(12,lb)=0
       endif
       elseif(jf.eq.3) then
       if(tneigh_ref_lev(1).ne.l_ref) mark_edge(1,lb)=0
       if(tneigh_ref_lev(2).ne.l_ref) mark_edge(3,lb)=0
       if(ndim.eq.3) then
       if(tneigh_ref_lev(5).ne.l_ref) mark_edge(5,lb)=0
       if(tneigh_ref_lev(6).ne.l_ref) mark_edge(7,lb)=0
       endif
       elseif(jf.eq.4) then
       if(tneigh_ref_lev(1).ne.l_ref) mark_edge(2,lb)=0
       if(tneigh_ref_lev(2).ne.l_ref) mark_edge(4,lb)=0
       if(ndim.eq.3) then
       if(tneigh_ref_lev(5).ne.l_ref) mark_edge(6,lb)=0
       if(tneigh_ref_lev(6).ne.l_ref) mark_edge(8,lb)=0
       endif
       elseif(jf.eq.5) then
       if(tneigh_ref_lev(1).ne.l_ref) mark_edge(9,lb)=0
       if(tneigh_ref_lev(2).ne.l_ref) mark_edge(11,lb)=0
       if(tneigh_ref_lev(3).ne.l_ref) mark_edge(5,lb)=0
       if(tneigh_ref_lev(4).ne.l_ref) mark_edge(6,lb)=0
       elseif(jf.eq.6) then
       if(tneigh_ref_lev(1).ne.l_ref) mark_edge(10,lb)=0
       if(tneigh_ref_lev(2).ne.l_ref) mark_edge(12,lb)=0
       if(tneigh_ref_lev(3).ne.l_ref) mark_edge(7,lb)=0
       if(tneigh_ref_lev(4).ne.l_ref) mark_edge(8,lb)=0
       endif

       endif
       enddo ! loop over faces


! Any edges on this block which are still marked need a diagonal patch.
! Note that in the shmem_gets below, we can always assume that a
! neighbor block exists, since the edge would not have been marked
! earlier if that was not so.

! Loop over the edges on this block.
       do ie=1,nbedges

       if(mark_edge(ie,lb).eq.1) then


! The edge data on the neighboring faces can be assumed to have been averaged
! correctly from the refined diagonal blocks.

! Now copy over the edge data from one of the neighbors.
       if(ie.eq.1) then
       jf=1
                      remote_pe    = neigh(2,jf,lb)
                      remote_block = neigh(1,jf,lb)


#if N_DIM >= 2
       do k=klo,kup
!                        call shmem_real_get
!     . (
!     .         receive(1,k),
!     .                    bedge_facex_z(1,2,1+nguard,
!     .           k,remote_block),
!     .                    nedgevar,remote_pe)
                        bedge_facex_z(:,1,1+nguard,k,lb)= & 
     &         receive(:,k)
       enddo
                        bedge_facey_z(:,1+nguard,1,klo:kup,lb)= & 
     &                    bedge_facex_z(:,1,1+nguard,klo:kup,lb)
#endif


       elseif(ie.eq.2) then
                        jf=1
                        remote_pe    = neigh(2,jf,lb)
                        remote_block = neigh(1,jf,lb)

#if N_DIM >= 2
       do k=klo,kup
!                        call shmem_real_get
!     . (
!     .                    receive(1,k),
!     .                    bedge_facex_z(1,2,1+nguard+nyb,
!     .                                  k,remote_block),
!     .                    nedgevar,remote_pe)
                        bedge_facex_z(:,1,1+nguard+nyb,k,lb)= & 
     &                    receive(:,k)
       enddo
                        bedge_facey_z(:,1+nguard,2,klo:kup,lb)= & 
     &                    bedge_facex_z(:,1,1+nguard+nyb,klo:kup,lb)
#endif

                elseif(ie.eq.3) then
                        jf=2
                        remote_pe    = neigh(2,jf,lb)
                        remote_block = neigh(1,jf,lb)

#if N_DIM >= 2
       do k=klo,kup
!                        call shmem_real_get
!     . (
!     .                    receive(1,k),
!     .                    bedge_facex_z(1,1,1+nguard,
!     .                                  k,remote_block),
!     .                    nedgevar,remote_pe)
                        bedge_facex_z(:,2,1+nguard,k,lb)= & 
     &                    receive(:,k)
       enddo
                        bedge_facey_z(:,1+nguard+nxb,1,klo:kup,lb)= & 
     &                    bedge_facex_z(:,2,1+nguard,klo:kup,lb)
#endif


                elseif(ie.eq.4) then
                        jf=2
                        remote_pe    = neigh(2,jf,lb)
                        remote_block = neigh(1,jf,lb)

#if N_DIM >= 2
       do k=klo,kup
!                        call shmem_real_get
!     . (
!     .                    receive(1,k),
!     .                    bedge_facex_z(1,1,1+nguard+nyb,
!     .                                  k,remote_block),
!     .                    nedgevar,remote_pe)
                        bedge_facex_z(:,2,1+nguard+nyb,k,lb)= & 
     &                    receive(:,k)
       enddo
                        bedge_facey_z(:,1+nguard+nxb,2,klo:kup,lb)= & 
     &                    bedge_facex_z(:,2,1+nguard+nyb,klo:kup,lb)
#endif


                elseif(ie.eq.5) then
                        jf=3
                        remote_pe    = neigh(2,jf,lb)
                        remote_block = neigh(1,jf,lb)

#if N_DIM >= 2
       do i=ilo,iup
!                        call shmem_real_get
!     . (
!     .                    receive(1,i),
!     .                    bedge_facey_x(1,i,2,
!     .                                  1+nguard*k3d,remote_block),
!     .                    nedgevar,remote_pe)
                        bedge_facey_x(:,i,1,klo,lb)= & 
     &                    receive(:,i)
       enddo
                        bedge_facez_x(:,ilo:iup,1+nguard*k3d,1,lb)= & 
     &                    bedge_facey_x(:,ilo:iup,1,klo,lb)
#endif


                elseif(ie.eq.6) then
                        jf=4
                        remote_pe    = neigh(2,jf,lb)
                        remote_block = neigh(1,jf,lb)

#if N_DIM >= 2
       do i=ilo,iup
!                        call shmem_real_get
!     . (
!     .                    receive(1,i),
!     .                    bedge_facey_x(1,i,1,
!     .                                  klo,remote_block),
!     .                    nedgevar,remote_pe)
                        bedge_facey_x(:,i,2,klo,lb)= & 
     &                    receive(:,i)
       enddo
                        bedge_facez_x(:,ilo:iup,1+nguard+nyb,1,lb)= & 
     &                    bedge_facey_x(:,ilo:iup,2,klo,lb)
#endif


                elseif(ie.eq.7) then
                        jf=3
                        remote_pe    = neigh(2,jf,lb)
                        remote_block = neigh(1,jf,lb)

#if N_DIM >= 2
       do i=ilo,iup
!                        call shmem_real_get
!     . (
!     .                    receive(1,i),
!     .                    bedge_facey_x(1,i,2,
!     .                                  kup,remote_block),
!     .                    nedgevar,remote_pe)
                        bedge_facey_x(:,i,1,kup,lb)= & 
     &                    receive(:,i)
       enddo
                        bedge_facez_x(:,ilo:iup,1+nguard,2,lb)= & 
     &                    bedge_facey_x(:,ilo:iup,1,kup,lb)
#endif


                elseif(ie.eq.8) then
                        jf=4
                        remote_pe    = neigh(2,jf,lb)
                        remote_block = neigh(1,jf,lb)

#if N_DIM >= 2
       do i=ilo,iup
!                        call shmem_real_get
!     . (
!     .                    receive(1,i),
!     .                    bedge_facey_x(1,i,1,
!     .                                  kup,remote_block),
!     .                    nedgevar,remote_pe)
                        bedge_facey_x(:,i,2,kup,lb)= & 
     &                    receive(:,i)
       enddo
                        bedge_facez_x(:,ilo:iup,1+nguard+nyb,2,lb)= & 
     &                    bedge_facey_x(:,ilo:iup,2,kup,lb)
#endif


                elseif(ie.eq.9) then
                        jf=1
                        remote_pe    = neigh(2,jf,lb)
                        remote_block = neigh(1,jf,lb)

#if N_DIM >= 2
       do j=jlo,jup
!                        call shmem_real_get
!     . (
!     .                    receive(1,j),
!     .                    bedge_facex_y(1,2,j,
!     .                                  klo,remote_block),
!     .                    nedgevar,remote_pe)
                        bedge_facex_y(:,1,j,klo,lb)= & 
     &                    receive(:,j)
       enddo
                        bedge_facez_y(:,1+nguard,jlo:jup,1,lb)= & 
     &                    bedge_facex_y(:,1,jlo:jup,klo,lb)
#endif


                elseif(ie.eq.10) then
                        jf=1
                        remote_pe    = neigh(2,jf,lb)
                        remote_block = neigh(1,jf,lb)

#if N_DIM >= 2
       do j=jlo,jup
!                        call shmem_real_get
!     . (
!     .                    receive(1,j),
!     .                    bedge_facex_y(1,2,j,
!     .                                  kup,remote_block),
!     .                    nedgevar,remote_pe)
                        bedge_facex_y(:,1,j,kup,lb)= & 
     &                    receive(:,j)
       enddo
                        bedge_facez_y(:,1+nguard,jlo:jup,2,lb)= & 
     &                    bedge_facex_y(:,1,jlo:jup,kup,lb)
#endif


                elseif(ie.eq.11) then
                        jf=2
                        remote_pe    = neigh(2,jf,lb)
                        remote_block = neigh(1,jf,lb)


#if N_DIM >= 2
       do j=jlo,jup
!                        call shmem_real_get
!     . (
!     .                    receive(1,j),
!     .                    bedge_facex_y(1,1,j,
!     .                                  klo,remote_block),
!     .                    nedgevar,remote_pe)
                        bedge_facex_y(:,2,j,klo,lb)= & 
     &                    receive(:,j)
       enddo
                        bedge_facez_y(:,1+nguard+nxb,jlo:jup,1,lb)= & 
     &                    bedge_facex_y(:,2,jlo:jup,klo,lb)
#endif


                elseif(ie.eq.12) then
                        jf=2
                        remote_pe    = neigh(2,jf,lb)
                        remote_block = neigh(1,jf,lb)

#if N_DIM >= 2
       do j=jlo,jup
!                        call shmem_real_get
!     . (
!     .                    receive(1,j),
!     .                    bedge_facex_y(1,1,j,
!     .                                  kup,remote_block),
!     .                    nedgevar,remote_pe)
                        bedge_facex_y(:,2,j,kup,lb)= & 
     &                    receive(:,j)
       enddo
                        bedge_facez_y(:,1+nguard+nxb,jlo:jup,2,lb)= & 
     &                    bedge_facex_y(:,2,jlo:jup,kup,lb)
#endif

       endif


       endif

       enddo ! loop over edges

      endif

      enddo
      endif


!        call shmem_barrier_all()


      return
      end
