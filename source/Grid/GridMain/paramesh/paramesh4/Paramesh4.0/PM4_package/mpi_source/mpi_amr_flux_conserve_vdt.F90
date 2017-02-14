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

!!****fi mpi_source/amr_flux_conserve_vdt
!! NAME
!!
!!   amr_flux_conserve_vdt
!!
!! SYNOPSIS
!!
!!   call amr_flux_conserve_vdt (mype, nsub)
!!
!!   call amr_flux_conserve_vdt (integer, integer)
!!
!! ARGUMENTS
!!   
!!   integer, intent(in) :: mype          
!!     The calling processor number.
!!
!!   integer, intent(in) :: nsub          
!!     The current time subcycle. If this is 1 then this info is used to 
!!     reset the temporary boundary flux arrays to 0.
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!
!! USES
!! 
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!
!! CALLS
!! 
!!   amr_restrict_bnd_data_vdt
!!   mpi_amr_comm_setup
!!    
!! RETURNS
!!
!!   Does not return anything.  Upon exit the fluxes stored in the arrays
!!   flux_x, flux_y, or flux_z are corrected at jumps in refinement.  This
!!   is either an averaging proceedure or a sum as selected by the user
!!   by adjusting the preprocessor variables which control this behaviour
!!   in 'paramesh_preprocessor.fh'.  This can also be controled through the 
!!   equivalent logical variables read at runtime from the file 
!!   'amr_runtime_parameters' if the 'LIBRARY' preprocessor flag is defined.
!!
!! DESCRIPTION
!!
!!   This routine gets block boundary data from neighbors who are
!!   parents of leaf blocks. This is required in flux conserving schemes
!!   where the coarser block needs to use the same fluxes and mean pressures
!!   as will be used on the finer blocks across their shared boundary.
!!
!!   The data structure used to store and pass this data is defined
!!   in the include file 'block_boundary_data.fh' which can be included
!!   in 'physicaldata.fh'.
!!
!!   This version is used when variable timesteps are allowed across the
!!   blocks in the computation.
!!
!! NOTES
!!  
!!   This routine is NOT user callable.  The user interacts with this routine
!!   by calling 'amr_flux_conserve' which, in turn, calls this routine.
!!
!!   THIS ROUTINE HAS NOT BEEN TESTED - July 14,2000
!!
!! AUTHORS
!!
!!   Peter MacNeice (1997) 
!!
!!***

      subroutine amr_flux_conserve_vdt(mype,nsub)


      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton
      use paramesh_comm_data

      use paramesh_interfaces, only : amr_restrict_bnd_data_vdt, & 
     &                                amr_mpi_find_blk_in_buffer
      use paramesh_mpi_interfaces, only : mpi_amr_comm_setup, & 
     &                                    mpi_set_message_limits

      implicit none

      include 'mpif.h'

      integer, intent(in)  ::  mype,nsub

!------------------------------------
! local variables

      integer :: remote_pe,remote_block
      integer :: remote_pe2,remote_block2
      integer,save ::  anodetype(1)
      integer :: cnodetype
      logical :: lfound
      save       cnodetype

      integer :: tag_offset,nprocs,ierr

      logical :: lcc ,lfc,lec,lnc
      logical :: lguard,lprolong,lflux,ledge,lrestrict,lfulltree
      integer :: iopt
      integer :: lb, lcycle, jf, iblk, iface
      integer :: i, j, k, ia0, ib0, ja0, jb0, ka0, kb0
      integer :: dtype, vtype, index, index0
      integer :: ia, ib, ja, jb, ka, kb
      integer :: n

      real :: phase0, phase1


!------------------------------------

      if (var_dt) then

      lcc = .false.
      lfc = .false.
      lec = .false.
      lnc = .false.
      iopt = 1


      if(lnblocks.gt.0) then
      do lb = 1,lnblocks

! Is this a parent of at least one leaf block ?
      if(nodetype(lb).eq.2) then

! Set timestep phases for the current block, and for the next finer level.
        lcycle = loc_cycle(lrefine(lb))
        phase0 = phase_dt(lrefine(lb))
        phase1 = phase_dt(lrefine(lb)+1)

! At start of the current blocks timestep zero out the arrays used to 
! accumulate boundary fluxes from its children.
        if(lcycle.eq.1) then
           ttflux_x(:,:,:,:,lb) = 0.
           if(ndim.ge.2) ttflux_y(:,:,:,:,lb) = 0.
           if(ndim.eq.3) ttflux_z(:,:,:,:,lb) = 0.
        endif

      endif
      enddo
      endif
!------------------------------------

      Call MPI_COMM_SIZE(amr_mpi_meshComm, nprocs, ierr)
      tag_offset = 100


! Note, both lflux and lrestrict are true so that the fluxes
! are acquired which are needed in the restriction operation.
      lguard    = .false.
      lprolong  = .false.
      lflux     = .true.
      ledge     = .false.
      lrestrict = .true.
      lrestrict = .false.
      call mpi_amr_comm_setup(mype,nprocs,lguard,lprolong, & 
     &                        lflux,ledge,lrestrict,lfulltree, & 
     &                        iopt,lcc,lfc,lec,lnc,tag_offset)


! Leaf blocks which have completed their timestep provide reduced 
! boundary data to their parents.
! Fluxes are accumulated in the ttflux_ arrays.
      call amr_restrict_bnd_data_vdt(mype)

!------------------------------------


! Parents who have completed their timestep and border a leaf block
! update their fluxes.
      do lb = 1,lnblocks


! Is this a parent block of at least one leaf node?
      if((nodetype(lb).eq.2).and.ldtcomplete(lb)) then

! If yes then cycle through its neighbors.
        do iface=1,nfaces

! If this neighbor is a leaf block or an external boundary then 
! replace fluxes with restricted fluxes.
          cnodetype = 1
          if(neigh(1,iface,lb).ge.1) then
            remote_pe    = neigh(2,iface,lb)
            remote_block = neigh(1,iface,lb)

! if (remote_block,remote_pe) is not a local block then it must have a
! local copy available in the buffer space at the end of the local
! block list.

          if(remote_pe.ne.mype) then

            do iblk = strt_buffer,last_buffer
              if(remote_block.eq.laddress(1,iblk).and. & 
     &             remote_pe .eq.laddress(2,iblk) ) then
                remote_block = iblk
                remote_pe    = mype
              endif
            enddo
          endif

            cnodetype = nodetype(remote_block)

          endif
          if(cnodetype.eq.1) then
            if(iface.eq.1) flux_x(:,1,:,:,lb)=ttflux_x(:,1,:,:,lb)
            if(iface.eq.2) flux_x(:,2,:,:,lb)=ttflux_x(:,2,:,:,lb)
            if(iface.eq.3) flux_y(:,:,1,:,lb)=ttflux_y(:,:,1,:,lb)
            if(iface.eq.4) flux_y(:,:,2,:,lb)=ttflux_y(:,:,2,:,lb)
            if(iface.eq.5) flux_z(:,:,:,1,lb)=ttflux_z(:,:,:,1,lb)
            if(iface.eq.6) flux_z(:,:,:,2,lb)=ttflux_z(:,:,:,2,lb)
          endif
        enddo
      endif
      enddo

      tag_offset = 100


      lguard    = .false.
      lprolong  = .false.
      lflux     = .true.
      ledge     = .false.
      lrestrict = .false.
      lfulltree = .false.
      call mpi_amr_comm_setup(mype,nprocs,lguard,lprolong, & 
     &                        lflux,ledge,lrestrict,lfulltree, & 
     &                        iopt,lcc,lfc,lec,lnc,tag_offset)

!------------------------------------

! cycle through the grid blocks on this processor
      if(lnblocks.gt.0) then
      do lb = 1,lnblocks

! Is this a leaf block and not at the original refinement level ?
      if(nodetype(lb).eq.1) then

! Has this block completed its timestep?
      if(ldtcomplete(lb)) then

! Cycle over the blocks faces
       do jf = 1,nfaces

          remote_pe = neigh(2,jf,lb)
          remote_block  = neigh(1,jf,lb)
          remote_pe2 = neigh(2,jf,lb)
          remote_block2  = neigh(1,jf,lb)
          cnodetype = 0
          lfound = .false.


          if(remote_block.gt.0) then

! (remote_block,remote_pe) may be a local block, a remote block,
! or it may not exist.
! If it is a local block then check its nodetype.
! If it is found in the list of remote blocks stored in buffer space
! then check its nodetype.
! If it is not found in either of these places, then set its nodetype
! to 0.
          if(remote_pe2.ne.mype) then

             lfound = .false.
             do iblk = strt_buffer,last_buffer
                if(remote_block2.eq.laddress(1,iblk).and. & 
     &               remote_pe2 .eq.laddress(2,iblk) ) then
                   remote_block2 = iblk
                   remote_pe2    = mype
                   lfound = .true.
                endif
             enddo
             
          elseif(remote_pe2.eq.mype) then

             lfound = .true.
             
          endif

! Is the neighbor to this face a parent of a leaf block?
          if(lfound) then
             cnodetype = nodetype(remote_block2)
          endif

          endif  !  end of remote_block if test

          if(cnodetype.eq.2) then

! If yes then copy the appropriate layer from its boundary variable data 

            if (remote_pe == mype .and. remote_block <= lnblocks) then

            if(jf.eq.1) then

               flux_x(1:nfluxes,1,:,:,lb) = flux_x(1:nfluxes,2,:,:,remote_block)

            elseif(jf.eq.2) then

               flux_x(1:nfluxes,2,:,:,lb) = flux_x(1:nfluxes,1,:,:,remote_block)

            elseif(jf.eq.3) then

               flux_y(1:nfluxes,:,1,:,lb) = flux_y(1:nfluxes,:,2,:,remote_block)

            elseif(jf.eq.4) then

               flux_y(1:nfluxes,:,2,:,lb) = flux_y(1:nfluxes,:,1,:,remote_block)

            elseif(jf.eq.5) then

               flux_z(1:nfluxes,:,:,1,lb) = flux_z(1:nfluxes,:,:,2,remote_block)

            elseif(jf.eq.6) then

               flux_z(1:nfluxes,:,:,2,lb) = flux_z(1:nfluxes,:,:,1,remote_block)

            endif      ! end of jf if


            else ! if (remote_pe


            call amr_mpi_find_blk_in_buffer(mype,remote_block, & 
     &              remote_pe,1,dtype,index0,lfound)
            vtype = 1
            call mpi_set_message_limits(dtype, & 
     &           ia0,ib0,ja0,jb0,ka0,kb0,vtype)

            index = index0 + 1

            if(dtype.eq.13.or.dtype.eq.15.or.dtype.eq.14) then

            ia = ia0
            ib = ib0
            ja = ja0
            jb = jb0
            ka = ka0
            kb = kb0

            if(dtype.eq.13) then
               ia = 1
               ib = 1
            elseif(dtype.eq.15) then
               ia = 2
               ib = 2
            elseif(dtype.eq.14) then
               ia = 1
               ib = 2
            endif

            if (jf == 1 .or. jf == 2) then

            do k = ka,kb
               do j = ja,jb
                  do i = ia,ib
                     do n=1,nfluxes
                        recvarxf(n,i,j,k) = temprecv_buf(index)
                        index  = index + 1
                     enddo
                  enddo
               enddo
            enddo

            end if

            end if

            if (ndim >= 2) then
               if(dtype.eq.11.or.dtype.eq.17.or.dtype.eq.14) then

                  ia = ia0
                  ib = ib0
                  ja = ja0
                  jb = jb0
                  ka = ka0
                  kb = kb0

                  if(dtype.eq.11) then
                     ja = 1
                     jb = 1
                  elseif(dtype.eq.17) then
                     ja = 2
                     jb = 2
                  elseif(dtype.eq.14) then
                     ja = 1
                     jb = 2
                  endif

                  if (jf == 3 .or. jf == 4) then

                  do k = ka,kb
                     do j = ja,jb
                        do i = ia,ib
                           do n=1,nfluxes
                              recvaryf(n,i,j,k) = & 
     &                             temprecv_buf(index)
                              index  = index + 1
                           enddo
                        enddo
                     enddo
                  enddo

                  endif

               endif
            endif

            if (ndim == 3) then
               if(dtype.eq.5.or.dtype.eq.23.or.dtype.eq.14) then

                  ia = ia0
                  ib = ib0
                  ja = ja0
                  jb = jb0
                  ka = ka0
                  kb = kb0

                  if(dtype.eq.5) then
                     ka = 1
                     kb = 1
                  elseif(dtype.eq.23) then
                     ka = 2
                     kb = 2
                  elseif(dtype.eq.14) then
                     ka = 1
                     kb = 2
                  endif

                  if (jf == 5 .or. jf == 6) then
                                       
                  do k = ka,kb
                     do j = ja,jb
                        do i = ia,ib
                           do n=1,nfluxes
                              recvarzf(n,i,j,k) = & 
     &                             temprecv_buf(index)
                              index  = index + 1
                           enddo
                        enddo
                     enddo
                  enddo

                  endif

               endif
            endif

            if(jf.eq.1) then

               flux_x(1:nfluxes,1,:,:,lb) = recvarxf(1:nfluxes,2,:,:)

            elseif(jf.eq.2) then

               flux_x(1:nfluxes,2,:,:,lb) = recvarxf(1:nfluxes,1,:,:)

            elseif(jf.eq.3) then

               flux_y(1:nfluxes,:,1,:,lb) = recvaryf(1:nfluxes,:,2,:)

            elseif(jf.eq.4) then

               flux_y(1:nfluxes,:,2,:,lb) = recvaryf(1:nfluxes,:,1,:)

            elseif(jf.eq.5) then

               flux_z(1:nfluxes,:,:,1,lb) = recvarzf(1:nfluxes,:,:,2)

            elseif(jf.eq.6) then

               flux_z(1:nfluxes,:,:,2,lb) = recvarzf(1:nfluxes,:,:,1)

            endif

            endif ! if (remote_pe

          endif        ! end of cnodetype if

       enddo

      endif                      ! end of ldtcomplete if test

      endif
      enddo
      endif

!------------------------------------

      endif

      return
      end subroutine amr_flux_conserve_vdt
