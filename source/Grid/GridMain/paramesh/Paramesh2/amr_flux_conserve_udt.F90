      subroutine amr_flux_conserve_udt(mype,idir)


! $RCSfile: amr_flux_conserve_udt.F90,v $
! $Revision: 1.2 $
! $Date: 2004/09/07 00:58:08 $


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
! This version is called when uniform timesteps are being used across
! the blocks in the computation.
!
!
! Written :     Peter MacNeice          February 1997
!------------------------------------------------------------------------

use physicaldata
      use tree
      implicit none
      include 'mpif.h'


!------------------------------------
! local variables

      integer mype,idir

      integer reqr(maxblocks)
      integer statr(MPI_STATUS_SIZE,maxblocks)
      integer remote_pe,remote_block
      integer cnodetype
      integer block_point(maxblocks)
      integer errorcode

      integer istart,iend,jstart,jend,kstart,kend
      integer jf,jsend,isg
      integer neighr,neighs

      integer len,ierr

      integer, parameter :: block_point_not_found = -99

      save      cnodetype

!------------------------------------

! all leaf blocks provide reduced boundary data to their parents
      call amr_restrict_bnd_data(mype,idir)

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

! cycle through the grid blocks on this processor
      if(lnblocks.gt.0) then

      block_point = block_point_not_found

! Cycle over the blocks faces
      do jf = istart,iend

        if (jf.eq.1) jsend = 2
        if (jf.eq.2) jsend = 1
        if (jf.eq.3) jsend = 4
        if (jf.eq.4) jsend = 3
        if (jf.eq.5) jsend = 6
        if (jf.eq.6) jsend = 5

        neighr = 0
        neighs = 0

        if (jf.eq.1.or.jf.eq.2) then

        len = len_block_bndx*nfluxes
        neighr = 0
        do isg = 1,lnblocks
          if (nodetype(isg).eq.1) then
            if (neigh_type(jf,isg).eq.2) then
              if(neigh(1,jf,isg).gt.-1) then
                if(neigh(2,jf,isg).ne.mype) then
                  neighr = neighr + 1
                  block_point(isg) = neighr + lnblocks
                  call MPI_IRECV(flux_x(1,1,1,1, & 
     &                 block_point(isg)), & 
     &                 len, & 
     &                 MPI_DOUBLE_PRECISION, & 
     &                 neigh(2,jf,isg), & 
     &                 neigh(1,jf,isg), & 
     &                 MPI_COMM_WORLD, & 
     &                 reqr(neighr), & 
     &                 ierr)
                end if
              end if
            end if
          end if
        end do
        
        if (lnblocks+neighr.gt.maxblocks) then
          print *,' ERROR memory overflow flux_conserve_udt ',mype
          call Driver_abortFlash("ERROR: memory overflow flux_conserve_udt")
        end if
        
! send messages if neighbor is off processor
      
        neighs = 0
        do isg = 1,lnblocks
          if (nodetype(isg).eq.2) then
            if (neigh_type(jsend,isg).eq.1) then
              if(neigh(1,jsend,isg).gt.-1) then
                if(neigh(2,jsend,isg).ne.mype) then
                  neighs = neighs + 1
                  call MPI_SSEND(flux_x(1,1,1,1, & 
     &                 isg), & 
     &                 len, & 
     &                 MPI_DOUBLE_PRECISION, & 
     &                 neigh(2,jsend,isg), &  ! PE TO SEND TO
     &                 isg,     &  ! THIS IS THE TAG
     &                 MPI_COMM_WORLD, & 
     &                 ierr)
                end if
              end if
            end if
          end if
        end do

        elseif (jf.eq.3.or.jf.eq.4) then

        len = len_block_bndy*nfluxes
        neighr = 0
        do isg = 1,lnblocks
          if (nodetype(isg).eq.1) then
            if (neigh_type(jf,isg).eq.2) then
              if(neigh(1,jf,isg).gt.-1) then
                if(neigh(2,jf,isg).ne.mype) then
                  neighr = neighr + 1
                  block_point(isg) = neighr + lnblocks
                  call MPI_IRECV(flux_y(1,1,1,1, & 
     &                 block_point(isg)), & 
     &                 len, & 
     &                 MPI_DOUBLE_PRECISION, & 
     &                 neigh(2,jf,isg), & 
     &                 neigh(1,jf,isg), & 
     &                 MPI_COMM_WORLD, & 
     &                 reqr(neighr), & 
     &             ierr)
                end if
              end if
            end if
          end if
        end do
        
        if (lnblocks+neighr.gt.maxblocks) then
          print *,' ERROR memory overflow in flux_conserve_ut ',mype
          call Driver_abortFlash("ERROR: memory overflow in flux_conserve_ut")
        end if
        
! send messages if neighbor is off processor
      
        neighs = 0
        do isg = 1,lnblocks
          if (nodetype(isg).eq.2) then
            if (neigh_type(jsend,isg).eq.1) then
              if(neigh(1,jsend,isg).gt.-1) then
                if(neigh(2,jsend,isg).ne.mype) then
                  neighs = neighs + 1
                  call MPI_SSEND(flux_y(1,1,1,1, & 
     &                 isg), & 
     &                 len, & 
     &                 MPI_DOUBLE_PRECISION, & 
     &                 neigh(2,jsend,isg), &  ! PE TO SEND TO
     &                 isg,     &  ! THIS IS THE TAG
     &                 MPI_COMM_WORLD, & 
     &                 ierr)
                end if
              end if
            end if
          end if
        end do

        elseif (jf.eq.5.or.jf.eq.6) then

        len = len_block_bndz*nfluxes
        neighr = 0
        do isg = 1,lnblocks
          if (nodetype(isg).eq.1) then
            if (neigh_type(jf,isg).eq.2) then
              if(neigh(1,jf,isg).gt.-1) then
                if(neigh(2,jf,isg).ne.mype) then
                  neighr = neighr + 1
                  block_point(isg) = neighr + lnblocks
                  call MPI_IRECV(flux_z(1,1,1,1, & 
     &                 block_point(isg)), & 
     &                 len, & 
     &                 MPI_DOUBLE_PRECISION, & 
     &                 neigh(2,jf,isg), & 
     &                 neigh(1,jf,isg), & 
     &                 MPI_COMM_WORLD, & 
     &                 reqr(neighr), & 
     &                 ierr)
                end if
              end if
            end if
          end if
        end do
        
        if (lnblocks+neighr.gt.maxblocks) then
          print *,' ERROR memory overflow in flux_conserve_ut ',mype
          call Driver_abortFlash("ERROR: memory overflow in flux_conserve_ut")
        end if
        
! send messages if neighbor is off processor
      
        neighs = 0
        do isg = 1,lnblocks
          if (nodetype(isg).eq.2) then
            if (neigh_type(jsend,isg).eq.1) then
              if(neigh(1,jsend,isg).gt.-1) then
                if(neigh(2,jsend,isg).ne.mype) then
                  neighs = neighs + 1
                  call MPI_SSEND(flux_z(1,1,1,1, & 
     &                 isg), & 
     &                 len, & 
     &                 MPI_DOUBLE_PRECISION, & 
     &                 neigh(2,jsend,isg), &  ! PE TO SEND TO
     &                 isg,     &  ! THIS IS THE TAG
     &                 MPI_COMM_WORLD, & 
     &                 ierr)
                end if
              end if
            end if
          end if
        end do

      end if
        
      if (neighr.gt.0) then
        call MPI_WAITALL (neighr, reqr, statr, ierr)
      end if

      do isg = 1,lnblocks

! Is this a leaf block?
      if(nodetype(isg).eq.1) then

          remote_pe = neigh(2,jf,isg)
          remote_block = neigh(1,jf,isg)
          if (remote_pe.ne.mype) remote_block = block_point(isg)

! Is the neighbor to this face a parent of a leaf block?

          cnodetype = neigh_type(jf,isg)

          if(cnodetype.eq.2) then

          if ( remote_block == block_point_not_found ) then
             print *, '[amr_flux_conserve_udt] ERROR: Attempt to access unidentified remote block.'
        call Driver_abortFlash("[amr_flux_conserve_udt] ERROR: Attempt to access unidentified remote block.")
          end if

! If yes then copy the appropriate layer from its boundary variable data 

            if(jf.eq.1) then
              flux_x(1:nfluxes,1,:,:,isg) =  & 
     &             flux_x(1:nfluxes,2,:,:,remote_block)
            elseif(jf.eq.2) then
              flux_x(1:nfluxes,2,:,:,isg) =  & 
     &             flux_x(1:nfluxes,1,:,:,remote_block)
            elseif(jf.eq.3) then
              flux_y(1:nfluxes,:,1,:,isg) =  & 
     &             flux_y(1:nfluxes,:,2,:,remote_block)
            elseif(jf.eq.4) then
              flux_y(1:nfluxes,:,2,:,isg) =  & 
     &             flux_y(1:nfluxes,:,1,:,remote_block)
            elseif(jf.eq.5) then
              flux_z(1:nfluxes,:,:,1,isg) =  & 
     &             flux_z(1:nfluxes,:,:,2,remote_block)
            elseif(jf.eq.6) then
              flux_z(1:nfluxes,:,:,2,isg) =  & 
     &             flux_z(1:nfluxes,:,:,1,remote_block)
            endif

          endif

         endif
      
       enddo

      enddo

      endif

      return
      end
 
