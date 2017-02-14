      subroutine amr_guardcell_cc_c_to_f(mype,iopt,nlayers,idir, & 
     &     block_point,child_n)


! $RCSfile: amr_guardcell_cc_c_to_f.F90,v $
! $Revision: 1.2 $
! $Date: 2004/09/07 00:58:08 $


!------------------------------------------------------------------------
! Rewritten to accomodate even or odd block sizes but only 2nd order schemes
!
! This routine manages the transfer of guard cell data to all
! leaf blocks from any neighbors which they have at a coarser 
! resolution.
!
! Written :     Peter MacNeice          January 1997
!------------------------------------------------------------------------
!
! Arguments:
!       mype            local processor number
!       iopt            a switch to control which data source is to be used
!                       iopt=1 will use 'unk'
!                       iopt=2 will use 'work'
!       nlayers         the number of guard cell layers at each boundary
!
!------------------------------------

use physicaldata
        use tree
        use workspace
        use paramesh_interfaces
        use Driver_interface, ONLY : Driver_abortFlash
        implicit none
        include 'mpif.h'





        integer mype,iopt,nlayers,idir
        integer istart,iend,jstart,jend,kstart,kend

        integer neighr,neighs,isg
        integer ichi,jf,jblock,jchild,ich
        integer ioff,joff,koff
        integer ierr

!------------------------------------
! local arrays
      real recv(nvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
      integer block_point(maxblocks)
      integer child_n(2,mchild,maxblocks)
      integer ichild(2,mchild)
      integer reqr(maxblocks),reqcs(maxblocks)
      integer reqcr(maxblocks)
      integer statr(MPI_STATUS_SIZE,maxblocks)
      integer errorcode

      integer remote_pe,remote_block

!------------------------------------

      if (idir.eq.1) then
        istart = 1
        iend = 2
      elseif (idir.eq.2) then
        istart = 3
        iend = 4
      elseif (idir.eq.3) then
        istart = 5
        iend = 6
      else
        istart = 1
        iend = nfaces
      end if

! children receive messages from their parents
      if(lnblocks.gt.0) then

        if (iopt.eq.1) then

        neighr = 0
        do isg = 1,lnblocks
          if (nodetype(isg).eq.1) then
            if(parent(1,isg).gt.-1) then
              if(parent(2,isg).ne.mype) then
                neighr = neighr + 1
                if (lnblocks+neighr.gt.maxblocks) then
                   print *,' ERROR: memory overflow in amr_guardcell_cc_c_to_f for UNK data'
                   call Driver_abortFlash("ERROR: memory overflow in amr_guardcell_cc_c_to_f")
                end if
                call MPI_IRECV(unk(1,1,1,1,lnblocks+neighr), & 
     &               len_block, & 
     &               MPI_DOUBLE_PRECISION, & 
     &               parent(2,isg), & 
     &               isg, & 
     &               MPI_COMM_WORLD, & 
     &               reqr(neighr), & 
     &               ierr)
                call MPI_IRECV(child_n(1,1,neighr), & 
     &               2*mchild, & 
     &               MPI_INTEGER, & 
     &               parent(2,isg), & 
     &               isg, & 
     &               MPI_COMM_WORLD, & 
     &               reqcr(neighr), & 
     &               ierr)
                block_point(isg) = lnblocks+neighr
              end if
            end if
          end if
        end do

! parents send messages to their children if off processor

        neighs = 0
        do isg = 1,lnblocks
          if (nodetype(isg).eq.2) then
            do ichi = 1,nchild
              if (child_type(ichi,isg).eq.1) then
                if(child(2,ichi,isg).ne.mype) then
                  neighs = neighs + 1
                  call MPI_SSEND(unk(1,1,1,1,isg), & 
     &                 len_block, & 
     &                 MPI_DOUBLE_PRECISION, & 
     &                 child(2,ichi,isg), &  ! PE TO SEND TO
     &                 child(1,ichi,isg),     &  ! THIS IS THE TAG
     &                 MPI_COMM_WORLD, & 
     &                 ierr)
                  call MPI_SSEND(child(1,1,isg), & 
     &                 2*mchild, & 
     &                 MPI_INTEGER, & 
     &                 child(2,ichi,isg), &  ! PE TO SEND TO
     &                 child(1,ichi,isg),     &  ! THIS IS THE TAG
     &                 MPI_COMM_WORLD, & 
     &                 ierr)
                end if
              end if
            end do
          end if
        end do


        if (neighr.gt.0) then
          call MPI_WAITALL (neighr, reqr, statr, ierr)
          call MPI_WAITALL (neighr, reqcr, statr, ierr)
        end if

        else ! iopt = 2

        neighr = 0
        do isg = 1,lnblocks
          if (nodetype(isg).eq.1) then
            if(parent(1,isg).gt.-1) then
              if(parent(2,isg).ne.mype) then
                neighr = neighr + 1
                if (lnblocks+neighr.gt.maxblocks) then
                   print *,' ERROR: memory overflow in amr_guardcell_cc_c_to_f for WORK data'
                   call Driver_abortFlash("ERROR: memory overflow in amr_guardcell_cc_c_to_f")
                end if
                call MPI_IRECV(work(1,1,1,lnblocks+neighr,1), & 
     &               len_wblock, & 
     &               MPI_DOUBLE_PRECISION, & 
     &               parent(2,isg), & 
     &               isg, & 
     &               MPI_COMM_WORLD, & 
     &               reqr(neighr), & 
     &               ierr)
                call MPI_IRECV(child_n(1,1,neighr), & 
     &               2*mchild, & 
     &               MPI_INTEGER, & 
     &               parent(2,isg), & 
     &               isg, & 
     &               MPI_COMM_WORLD, & 
     &               reqcr(neighr), & 
     &               ierr)
                block_point(isg) = lnblocks+neighr
              end if
            end if
          end if
        end do

! parents send messages to their children if off processor

        neighs = 0
        do isg = 1,lnblocks
          if (nodetype(isg).eq.2) then
            do ichi = 1,nchild
              if (child_type(ichi,isg).eq.1) then
                if(child(2,ichi,isg).ne.mype) then
                  neighs = neighs + 1
                  call MPI_SSEND(work(1,1,1,isg,1), & 
     &                 len_wblock, & 
     &                 MPI_DOUBLE_PRECISION, & 
     &                 child(2,ichi,isg), &  ! PE TO SEND TO
     &                 child(1,ichi,isg),     &  ! THIS IS THE TAG
     &                 MPI_COMM_WORLD, & 
     &                 ierr)
                  call MPI_SSEND(child(1,1,isg), & 
     &                 2*mchild, & 
     &                 MPI_INTEGER, & 
     &                 child(2,ichi,isg), &  ! PE TO SEND TO
     &                 child(1,ichi,isg),     &  ! THIS IS THE TAG
     &                 MPI_COMM_WORLD, & 
     &                 ierr)
                end if
              end if
            end do
          end if
        end do

        if (neighr.gt.0) then
          call MPI_WAITALL (neighr, reqr, statr, ierr)
          call MPI_WAITALL (neighr, reqcr, statr, ierr)
        end if

        end if

      end if

! cycle through the grid blocks on this processor

      if(lnblocks.gt.0) then
      do isg = 1,lnblocks

! Is this a leaf block ?
      if(nodetype(isg).eq.1) then

! Cycle over the blocks faces
       do jf = istart,iend

! Is the neighbor to this face at coarser resolution?
       if(neigh(1,jf,isg).eq.-1) then

! If yes then copy the appropriate layers from its parents data 
! into a temporary buffer array.

! Does the parent reside on this processor?
       if(parent(2,isg).eq.mype) then

! copy required layers from parents data - note if we use the standard
! numerical recipes prolongation operator below, then we must copy 3 layers
! here, ie the guardcell layers and 2 immediate interior layers.
         jblock = parent(1,isg)
         
         if(iopt.eq.2) then
           recv1(:,:,:) = work(:,:,:,jblock,1)
         endif

         ichild(:,:)=child(:,:,jblock)

       else

         remote_pe = parent(2,isg)
         remote_block  = parent(1,isg)

         if(iopt.eq.1) then
           jblock = block_point(isg)
         elseif(iopt.eq.2) then
           recv1(:,:,:) = work(:,:,:,block_point(isg),1)
         endif

         ichild(:,:) = child_n(:,:,block_point(isg)-lnblocks)

       endif

! identify which child of its parent this leaf block represents
       jchild = 0
       do ich=1,nchild
         if( (ichild(1,ich).eq.isg).and. & 
     &        (ichild(2,ich).eq.mype) ) jchild=ich
       enddo
       
! compute the offset in the parent block appropriate for this child
       ioff = mod(jchild-1,2)*nxb/2
       joff = mod((jchild-1)/2,2)*nyb/2
       koff = mod((jchild-1)/4,2)*nzb/2

! interpolate(prolongate) data from the parent to the child
       if(iopt.eq.1) then
         call amr_prolong_unk_fun(unk(:,:,:,:,jblock), & 
     &        isg,ioff,joff, & 
     &        koff,jf,mype, jblock)
       elseif(iopt.eq.2) then
         call amr_prolong_work_fun(ilw,iuw,jlw,juw, & 
     &        klw,kuw,nlayers,isg,ioff,joff,koff, & 
     &        jf,mype)
       endif
       
      
      endif
      
      enddo

      endif
      enddo
      endif
      
      return
      end
