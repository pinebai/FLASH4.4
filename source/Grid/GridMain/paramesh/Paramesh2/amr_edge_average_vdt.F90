      subroutine amr_edge_average_vdt(mype,nsub)


! $RCSfile: amr_edge_average_vdt.F90,v $
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
! Arguments:
!       mype            processor number
!       nsub            current time subcycle. If this is 1 then this
!                       info is used to reset the temporary boundary edge
!                       data arrays to 0.
!
!
! Written :     Peter MacNeice          August 1997
!------------------------------------------------------------------------

use physicaldata
      use tree
      implicit none


      integer mype,nsub

!------------------------------------
! local variables

      integer ng_off
      parameter(ng_off=nguard+iface_off)

        integer remote_pe,remote_block
        integer cnodetype,cchild(2),ich
        logical lnodetime
        save    lnodetime,cchild,cnodetype

!------------------------------------

#ifdef VAR_DT


        if(lnblocks.gt.0) then
        do isg = 1,lnblocks

! Is this a leaf block ?
        if(nodetype(isg).eq.1) then

! At start of overall timestep zero out the temporary arrays used to store
! boundary fluxes.
        if(nsub.eq.1) then
                tedge_facex_y(:,:,:,:,isg) = 0.
                tedge_facey_x(:,:,:,:,isg) = 0.
                if(ndim.eq.3) then
                tedge_facex_z(:,:,:,:,isg) = 0.
                tedge_facey_z(:,:,:,:,isg) = 0.
                tedge_facez_x(:,:,:,:,isg) = 0.
                tedge_facez_y(:,:,:,:,isg) = 0.
                endif
        endif

        endif
        enddo
        endif
!      call shmem_barrier_all()

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
!                if(remote_block.gt.0) call shmem_get(cnodetype,
!     .       nodetype(remote_block),1,remote_pe)
       if(cnodetype.eq.2) then
! If so, locate one of its children facing the current block
! (ie note that if jf is odd you can always use child nchild, while if
! jf is even you can always use child 1 - these choices work in both
! 2D and 3D) and find its timestep state.
                        ich = 1+(nchild-1)*mod(jf,2)
!                        call shmem_get(cchild,
!     .                   child(1,ich,remote_block),2,remote_pe)
!                        call shmem_get(lnodetime,
!     .                          ldtcomplete(cchild(1)),1,cchild(2))


! If yes then copy the appropriate layer from its boundary variable data 

       if(jf.eq.1) then

!       call shmem_get(recvarx1e(1,1,1,1),
!     .               bedge_facex_y(1,1,1,1,remote_block),
!     .       len_block_ex*nedges,remote_pe)
                        if(lnodetime) then
                                tbedge_facex_y(:,2,:,:,isg)= & 
     &                                  tbedge_facex_y(:,2,:,:,isg)+ & 
     &                                          recvarx1e(:,2,:,:)
                                bedge_facex_y(:,1,:,:,isg)= & 
     &                                  tbedge_facex_y(:,2,:,:,isg)
                        else
       bedge_facex_y(:,1,:,:,isg)= & 
     &       recvarx1e(:,2,:,:)
                        endif
                        if(ldtcomplete(isg)) & 
     &                          tbedge_facex_y(:,2,:,:,isg)= 0.

       if(ndim.eq.3) then
!       call shmem_get(recvarx1e(1,1,1,1),
!     .               bedge_facex_z(1,1,1,1,remote_block),
!     .       len_block_ex*nedges,remote_pe)
                        if(lnodetime) then
                                tbedge_facex_z(:,2,:,:,isg)= & 
     &                                  tbedge_facex_z(:,2,:,:,isg)+ & 
     &                                          recvarx1e(:,2,:,:)
                                bedge_facex_z(:,1,:,:,isg)= & 
     &                                  tbedge_facex_z(:,2,:,:,isg)
                        else
       bedge_facex_z(:,1,:,:,isg)= & 
     &       recvarx1e(:,2,:,:)
       endif
                        if(ldtcomplete(isg)) & 
     &                          tbedge_facex_z(:,2,:,:,isg)= 0.
       endif



       elseif(jf.eq.2) then

!       call shmem_get(recvarx1e(1,1,1,1),
!     .               bedge_facex_y(1,1,1,1,remote_block), & 
!     .       len_block_ex*nedges,remote_pe)
                        if(lnodetime) then
                                tbedge_facex_y(:,1,:,:,isg)= & 
     &                                  tbedge_facex_y(:,1,:,:,isg)+ & 
     &                                          recvarx1e(:,1,:,:)
                                bedge_facex_y(:,2,:,:,isg)= & 
     &                                  tbedge_facex_y(:,1,:,:,isg)
                        else
       bedge_facex_y(:,2,:,:,isg)= & 
     &       recvarx1e(:,1,:,:)
                        endif
                        if(ldtcomplete(isg)) & 
     &                          tbedge_facex_y(:,1,:,:,isg)= 0.



       if(ndim.eq.3) then
!       call shmem_get(recvarx1e(1,1,1,1),
!     .               bedge_facex_z(1,1,1,1,remote_block),
!     .       len_block_ex*nedges,remote_pe)
                        if(lnodetime) then
                                tbedge_facex_z(:,1,:,:,isg)= & 
     &                                  tbedge_facex_z(:,1,:,:,isg)+ & 
     &                                          recvarx1e(:,1,:,:)             
                                bedge_facex_z(:,2,:,:,isg)= & 
     &                                  tbedge_facex_z(:,1,:,:,isg)
                        else
       bedge_facex_z(:,2,:,:,isg)= & 
     &       recvarx1e(:,1,:,:)
                        endif
                        if(ldtcomplete(isg)) & 
     &                          tbedge_facex_z(:,1,:,:,isg)= 0.

       endif

       elseif(jf.eq.3) then

       if(ndim.eq.3) then
!       call shmem_get(recvary1e(1,1,1,1),
!     .               bedge_facey_z(1,1,1,1,remote_block),
!     .       len_block_ey*nedges,remote_pe)
                        if(lnodetime) then
                                tbedge_facey_z(:,:,2,:,isg)= & 
     &                                  tbedge_facey_z(:,:,2,:,isg)+ & 
     &                                          recvary1e(:,:,2,:)
                                bedge_facey_z(:,:,1,:,isg)= & 
     &                                  tbedge_facey_z(:,:,2,:,isg)
                        else
       bedge_facey_z(:,:,1,:,isg)= & 
     &       recvary1e(:,:,2,:)
       endif
                        if(ldtcomplete(isg)) & 
     &                          tbedge_facey_z(:,:,2,:,isg)= 0.
                        endif


!       call shmem_get(recvary1e(1,1,1,1),
!     .               bedge_facey_x(1,1,1,1,remote_block),
!     .       len_block_ey*nedges,remote_pe)
                        if(lnodetime) then
                                tbedge_facey_x(:,:,2,:,isg)= & 
     &                                  tbedge_facey_x(:,:,2,:,isg)+ & 
     &                                          recvary1e(:,:,2,:)
                                bedge_facey_x(:,:,1,:,isg)= & 
     &                                  tbedge_facey_x(:,:,2,:,isg)
                        else
       bedge_facey_x(:,:,1,:,isg)= & 
     &       recvary1e(:,:,2,:)
                        endif
                        if(ldtcomplete(isg)) & 
     &                          tbedge_facey_x(:,:,2,:,isg)= 0.


       elseif(jf.eq.4) then

       if(ndim.eq.3) then
!       call shmem_get(recvary1e(1,1,1,1),
!     .               bedge_facey_z(1,1,1,1,remote_block),
!     .       len_block_ey*nedges,remote_pe)
                        if(lnodetime) then
                                tbedge_facey_z(:,:,1,:,isg)= & 
     &                                  tbedge_facey_z(:,:,1,:,isg)+ & 
     &                                          recvary1e(:,:,1,:)
                                bedge_facey_z(:,:,2,:,isg)= & 
     &                                  tbedge_facey_z(:,:,1,:,isg)
                        else
       bedge_facey_z(:,:,2,:,isg)= & 
     &       recvary1e(:,:,1,:)
                        endif
                        if(ldtcomplete(isg)) & 
     &                          tbedge_facey_z(:,:,1,:,isg)= 0.

       endif

!       call shmem_get(recvary1e(1,1,1,1),
!     .               bedge_facey_x(1,1,1,1,remote_block),
!     .       len_block_ey*nedges,remote_pe)
                        if(lnodetime) then
                                tbedge_facey_x(:,:,1,:,isg)= & 
     &                                  tbedge_facey_x(:,:,1,:,isg)+ & 
     &                                          recvary1e(:,:,1,:)
                                bedge_facey_x(:,:,2,:,isg)= & 
     &                                  tbedge_facey_x(:,:,1,:,isg)
                        else
       bedge_facey_x(:,:,2,:,isg)= & 
     &       recvary1e(:,:,1,:)
                        endif
                        if(ldtcomplete(isg)) & 
     &                          tbedge_facey_x(:,:,1,:,isg)= 0.



       elseif(jf.eq.5) then

!       call shmem_get(recvarz1e(1,1,1,1),
!     .               bedge_facez_x(1,1,1,1,remote_block),
!     .       len_block_ez*nedges,remote_pe)
                        if(lnodetime) then
                                tbedge_facez_x(:,:,:,2,isg)= & 
     &                                  tbedge_facez_x(:,:,:,2,isg)+ & 
     &                                          recvarz1e(:,:,:,2)
                                bedge_facez_x(:,:,:,1,isg)= & 
     &                                  tbedge_facez_x(:,:,:,2,isg)
                        else
       bedge_facez_x(:,:,:,1,isg)= & 
     &       recvarz1e(:,:,:,2)
                        endif
                        if(ldtcomplete(isg)) & 
     &                          tbedge_facez_x(:,:,:,2,isg)= 0.


!       call shmem_get(recvarz1e(1,1,1,1),
!     .               bedge_facez_y(1,1,1,1,remote_block),
!     .       len_block_ez*nedges,remote_pe)
                        if(lnodetime) then
                                tbedge_facez_y(:,:,:,2,isg)= & 
     &                                  tbedge_facez_y(:,:,:,2,isg)+ & 
     &                                          recvarz1e(:,:,:,2)
                                bedge_facez_y(:,:,:,1,isg)= & 
     &                                  tbedge_facez_y(:,:,:,2,isg)
                        else
       bedge_facez_y(:,:,:,1,isg)= & 
     &       recvarz1e(:,:,:,2) 
                        endif
                        if(ldtcomplete(isg)) & 
     &                          tbedge_facez_y(:,:,:,2,isg)= 0.



       elseif(jf.eq.6) then 

!       call shmem_get(recvarz1e(1,1,1,1),
!     .               bedge_facez_x(1,1,1,1,remote_block),
!     .       len_block_ez*nedges,remote_pe)
                        if(lnodetime) then
                                tbedge_facez_x(:,:,:,1,isg)= & 
     &                                  tbedge_facez_x(:,:,:,1,isg)+ & 
     &                                          recvarz1e(:,:,:,1)
                                bedge_facez_x(:,:,:,2,isg)= & 
     &                                  tbedge_facez_x(:,:,:,1,isg)
                        else
       bedge_facez_x(:,:,:,2,isg)= & 
     &       recvarz1e(:,:,:,1)
                        endif
                        if(ldtcomplete(isg)) & 
     &                          tbedge_facez_x(:,:,:,1,isg)= 0.


!       call shmem_get(recvarz1e(1,1,1,1),
!     .               bedge_facez_y(1,1,1,1,remote_block),
!     .       len_block_ez*nedges,remote_pe)
                        if(lnodetime) then
                                tbedge_facez_y(:,:,:,1,isg)= & 
     &                                  tbedge_facez_y(:,:,:,1,isg)+ & 
     &                                          recvarz1e(:,:,:,1)
                                bedge_facez_y(:,:,:,2,isg)= & 
     &                                  tbedge_facez_y(:,:,:,1,isg)
                        else
       bedge_facez_y(:,:,:,2,isg)= & 
     &       recvarz1e(:,:,:,1)
                        endif
                        if(ldtcomplete(isg)) & 
     &                          tbedge_facez_y(:,:,:,1,isg)= 0.



       endif


       endif

       enddo

      endif
      enddo
      endif

!      call shmem_barrier_all()

#endif

      return
      end
