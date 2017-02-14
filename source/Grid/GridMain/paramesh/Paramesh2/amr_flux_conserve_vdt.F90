      subroutine amr_flux_conserve_vdt(mype,nsub)


! $RCSfile: amr_flux_conserve_vdt.F90,v $
! $Revision: 1.1 $
! $Date: 2004/05/20 01:03:44 $


!------------------------------------------------------------------------
!
! This routine gets block boundary data from neighbors who are
! parents of leaf blocks. This is required in flux conserving schemes
! where the coarser block needs to use the same fluxes and mean pressures
! as will be used on the finer blocks across their shared boundary.
!
! The data structure used to store and pass this data is defined
! in the include file 'block_boundary_data.h' which can be included
! in 'physicaldata.h'.
!
! This version is used when variable timesteps are allowed across the
! blocks in the computation.
!
! Arguments:
!      mype          processor number
!      nsub          current time subcycle. If this is 1 then this
!                     info is used to reset the temporary boundary flux
!                     arrays to 0.
!
! Written :     Peter MacNeice          February 1997
!------------------------------------------------------------------------

use physicaldata
      use tree
      implicit none




!------------------------------------
! local variables

      integer mype,nsub

      integer remote_pe,remote_block
      integer cnodetype,cchild(2),ich
      logical lnodetime

      integer istart,iend,jstart,jend,kstart,kend
      integer jf,jsend,isg
      integer neighr,neighs
      integer len,ierr

      save    lnodetime,cchild,cnodetype

!------------------------------------


      if(lnblocks.gt.0) then
      do isg = 1,lnblocks

! Is this a leaf block ?
      if(nodetype(isg).eq.1) then

! At start of overall timestep zero out the temporary arrays used to store
! boundary fluxes.
        if(nsub.eq.1) then
           tflux_x(:,:,:,:,isg) = 0.
           if(ndim.ge.2) tflux_y(:,:,:,:,isg) = 0.
           if(ndim.eq.3) tflux_z(:,:,:,:,isg) = 0.
        endif

      endif
      enddo
      endif
!      call shmem_barrier_all()

! all leaf blocks provide reduced boundary data to their parents
      call amr_restrict_bnd_data(mype)


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
          if(remote_block.gt.0) then
!             call shmem_integer_get(cnodetype,
!     .                       nodetype(remote_block),1,remote_pe)
          endif


          if(cnodetype.eq.2) then
! If so, locate one of its children facing the current block
! (ie note that if jf is odd you can always use child nchild, while if
! jf is even you can always use child 1 - these choices work in both
! 2D and 3D) and find its timestep state.
             ich = 1+(nchild-1)*mod(jf,2)
!             call shmem_integer_get(cchild,
!     .                       child(1,ich,remote_block),2,remote_pe)
!             call shmem_logical_get(lnodetime,
!     .                          ldtcomplete(cchild(1)),1,cchild(2))


! If yes then copy the appropriate layer from its boundary variable data.

! In the following code blocks, after we get the appropriate data from
! remote_block, we either use it directly if the remote_block has not
! yet reached the end of its timestep, or if it has we accumulate it
! into the temporary storage tflux_(x/y/z). If the local block is finishing
! its timestep we finish by setting tflux_(x/y/z) to 0.

             if(jf.eq.1) then
!               call shmem_real_get
!     . (recvarx1(1,1,1,1),
!     .               flux_x(1,1,1,1,remote_block),
!     .               len_block_bndx*nfluxes,remote_pe)
               if(lnodetime) then
                 tflux_x(1:nfluxes,2,:,:,isg)= & 
     &                           tflux_x(1:nfluxes,2,:,:,isg)+ & 
     &                           recvarx1(1:nfluxes,2,:,:)
                 flux_x(1:nfluxes,1,:,:,isg)= & 
     &                           tflux_x(1:nfluxes,2,:,:,isg)
               else
                 flux_x(1:nfluxes,1,:,:,isg) =  & 
     &                           recvarx1(1:nfluxes,2,:,:)
               endif
               if(ldtcomplete(isg)) tflux_x(1:nfluxes,2,:,:,isg)= 0.

             elseif(jf.eq.2) then
!               call shmem_real_get
!     . (recvarx1(1,1,1,1),
!     .               flux_x(1,1,1,1,remote_block),
!     .               len_block_bndx*nfluxes,remote_pe)
               if(lnodetime) then
                 tflux_x(1:nfluxes,1,:,:,isg)= & 
     &                             tflux_x(1:nfluxes,1,:,:,isg)+ & 
     &                             recvarx1(1:nfluxes,1,:,:)
                 flux_x(1:nfluxes,2,:,:,isg)= & 
     &                             tflux_x(1:nfluxes,1,:,:,isg)
               else
                 flux_x(1:nfluxes,2,:,:,isg) =  & 
     &                             recvarx1(1:nfluxes,1,:,:)
               endif
               if(ldtcomplete(isg)) tflux_x(1:nfluxes,1,:,:,isg)= 0.

             elseif(jf.eq.3) then
!               call shmem_real_get
!     . (recvary1(1,1,1,1),
!     .               flux_y(1,1,1,1,remote_block),
!     .               len_block_bndy*nfluxes,remote_pe)
               if(lnodetime) then
                 tflux_y(1:nfluxes,:,2,:,isg)= & 
     &                             tflux_y(1:nfluxes,:,2,:,isg)+ & 
     &                             recvary1(1:nfluxes,:,2,:)
                 flux_y(1:nfluxes,:,1,:,isg) =  & 
     &                             tflux_y(1:nfluxes,:,2,:,isg)
               else
                 flux_y(1:nfluxes,:,1,:,isg) =  & 
     &                             recvary1(1:nfluxes,:,2,:)
               endif
               if(ldtcomplete(isg)) tflux_y(1:nfluxes,:,2,:,isg)= 0.

             elseif(jf.eq.4) then
!               call shmem_real_get
!     . (recvary1(1,1,1,1),
!     .               flux_y(1,1,1,1,remote_block),
!     .               len_block_bndy*nfluxes,remote_pe)
               if(lnodetime) then
                 tflux_y(1:nfluxes,:,1,:,isg)= & 
     &                             tflux_y(1:nfluxes,:,1,:,isg)+ & 
     &                             recvary1(1:nfluxes,:,1,:)
                 flux_y(1:nfluxes,:,2,:,isg)= & 
     &                             tflux_y(1:nfluxes,:,1,:,isg)
               else
                 flux_y(1:nfluxes,:,2,:,isg) =  & 
     &                             recvary1(1:nfluxes,:,1,:)
               endif
               if(ldtcomplete(isg)) tflux_y(1:nfluxes,:,1,:,isg)= 0.

             elseif(jf.eq.5) then
!               call shmem_real_get
!     . (recvarz1(1,1,1,1),
!     .               flux_z(1,1,1,1,remote_block),
!     .               len_block_bndz*nfluxes,remote_pe)
               if(lnodetime) then
                 tflux_z(1:nfluxes,:,:,2,isg)= & 
     &                             tflux_z(1:nfluxes,:,:,2,isg)+ & 
     &                             recvarz1(1:nfluxes,:,:,2)
                 flux_z(1:nfluxes,:,:,1,isg)= & 
     &                             tflux_z(1:nfluxes,:,:,2,isg)
               else
                 flux_z(1:nfluxes,:,:,1,isg)= & 
     &                             recvarz1(1:nfluxes,:,:,2)
               endif
               if(ldtcomplete(isg)) tflux_z(1:nfluxes,:,:,2,isg)= 0.

             elseif(jf.eq.6) then
!               call shmem_real_get
!     . (recvarz1(1,1,1,1),
!     .               flux_z(1,1,1,1,remote_block),
!     .               len_block_bndz*nfluxes,remote_pe)
               if(lnodetime) then
                 tflux_z(1:nfluxes,:,:,1,isg)= & 
     &                             tflux_z(1:nfluxes,:,:,1,isg)+ & 
     &                             recvarz1(1:nfluxes,:,:,1)
                 flux_z(1:nfluxes,:,:,2,isg)= & 
     &                             tflux_z(1:nfluxes,:,:,1,isg)
               else
                 flux_z(1:nfluxes,:,:,2,isg) =  & 
     &                             recvarz1(1:nfluxes,:,:,1)
               endif
               if(ldtcomplete(isg)) tflux_z(1:nfluxes,:,:,1,isg)= 0.

             endif


          endif

        enddo

      endif
      enddo
      endif

!      call shmem_barrier_all()


      return
      end
