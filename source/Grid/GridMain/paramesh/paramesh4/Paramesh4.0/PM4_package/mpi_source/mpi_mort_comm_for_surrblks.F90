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

!#define DEBUG
!#define DEBUGX

      subroutine mpi_mort_comm_for_surrblks & 
     &                   (mype,nprocs,tag_offset)



!
! DESIGN ISSUES :
!  
!  Some neighbors of parents are requested unnecessarily because we
!  cannot verify that the corresponding child neighbors actually
!  exist until after we have received morton lists.
!  Can we improve list of requested blocks by identifying neighbors
!  of parents which may not be needed?
!
!
!------------------------------------------------------------------------
!
! This routine calculates the morton number for each block on mype.
! It stores the result along with the refinement level of each block into
! the array mortonbnd, and distributes this array among all processors.
!
!
! Written :     Peter MacNeice  and Michael Gehmeyr          February 2000
!------------------------------------------------------------------------
!
! Arguments:
!      mype           rank of local processor
!
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton
      use paramesh_comm_data
      use paramesh_mpi_interfaces, only : compress_list, & 
     &                                    morton_neighbors

      implicit none

      include 'mpif.h'


      integer, intent(in)    ::  mype,nprocs
      integer, intent(inout) ::  tag_offset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables

      real    :: pbsize(3),pcoord(3)
      real    :: xmin, ymin, zmin
      real    :: xmax, ymax, zmax

      integer :: lb,i,j
      integer :: morton(6),level,jstack
      integer ::  lbfirst,lblast
      integer :: mort_neigh(6,3,3,3)
      integer :: pmort_neigh(6,3,3,3)
      real    :: pbndbox(2,3)
!     integer :: neigh_morts(6,3,npts_neigh),indx(npts_neigh)
      integer,dimension (:,:,:),allocatable:: neigh_morts
      integer,dimension (:,:,:),allocatable:: tneigh_morts
      integer,dimension (:)    ,allocatable:: indx
      integer :: istart,iend
      integer :: i_pe,j_pe,rem_block,rem_pe
      integer :: no_of_comm_procs
      integer :: ierrorcode,ierr,allocation_status,ierror
      integer :: interp_max_orderf,interp_max_ordere
      integer :: interp_max_order ,interp_max_ordern,interx
      integer :: no_of_remote_neighs
      integer :: max_no_to_be_received
      integer :: max_no_of_blocks
      integer :: no_of_comms_to_send
      integer :: istack, ioff, joff, koff, k, itag, ll, kk
      integer :: itemp, kstack, iprocs, isize, isrc, idest
      integer,dimension (:),  allocatable :: recvrequest
      integer,dimension (:,:),allocatable :: recvstatus
      integer :: nguarda 
      integer :: ii, jj

      logical :: lremote,lswap,lfound
      logical :: is_remote,is_found
      logical :: morton_greater_than
      logical :: morton_equal
      logical :: morton_less_than
      logical,save :: l_on_pe = .false.

      integer :: npts_neigh1,npts_neigh2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef DEBUG_FLOW_TRACE
      write(*,*) 'pe ',mype,' entered mpi_mort_comm_for_surrblks'
#endif /* DEBUG_FLOW_TRACE */

       nguarda = max(nguard,nguard_work)

       lbfirst = 1
       lblast  = lnblocks

       npts_neigh1 = npts_neigh
       npts_neigh2 = npts_neigh+100
       allocate(neigh_morts(6,3,npts_neigh2))

  
!
!
! This routine assumes that the grid blocks are ordered by morton
! number and that any blocks with different refinement levels but
! the same morton number are ordered from coarse to fine.

! mark morton data out of date.
!      morton_limits_set = .false.

! Find highest order interpolant used in prolongation, if using face, edge
! or corner data. If Muscl is used then this will be interpreted as order 1.
       interp_max_orderf = 0
       interp_max_ordere = 0
       interp_max_ordern = 0
       do i = 1,nfacevar
          interx = max(interp_mask_facex(i), & 
     &                 interp_mask_facey(i), & 
     &                 interp_mask_facez(i))
          if(interx.eq.20) interx = 1
          interp_max_orderf = max(interp_max_orderf,interx)
       enddo
       do i = 1,nvaredge
          interx = interp_mask_ec(i)
          if(interx.eq.20) interx = 1
          interp_max_ordere = max(interp_max_ordere,interx)
       enddo
       do i = 1,nvarcorn
          interx = interp_mask_nc(i)
          if(interx.eq.20) interx = 1
          interp_max_ordern = max(interp_max_ordern,interx)
       enddo
       interp_max_order = max(interp_max_orderf,interp_max_ordere, & 
     &                        interp_max_ordern)


!--------------------------------------------------
!

! Compute xmin,ymin,zmin,xmax,ymax,zmax or get them from storage
      xmin = grid_xmin
      ymin = grid_ymin
      zmin = grid_zmin
      xmax = grid_xmax
      ymax = grid_ymax
      zmax = grid_zmax

!
! Initializations
      no_of_comm_procs = 0
      no_of_remote_neighs = 0
      max_no_to_be_received = 0
      max_no_to_send = 0
      commatrix_send = 0
      commatrix_recv = 0
      pe_source = -1
      pe_destination = -1
!     neigh_morts = -1
      no_of_comms_to_send = 0

!--------------------------------------------------
!
! Step 3.
! Construct a list of potential neighbors of all blocks on this
! processor, and potential neighbors of their parents.
! Exclude any which are on this processor.

      istack = 0

#ifdef DEBUG
      write(*,*) 'xmin,ymin,zmin,xmax,ymax,zmax ', & 
     & xmin,ymin,zmin,xmax,ymax,zmax
#endif /* DEBUG */


!     do lb=1,lnblocks
      do lb=lbfirst,lblast

!-------------

! First get the possible neighbors of the current block
      mort_neigh = -1
      pmort_neigh = -1
      call morton_neighbors(xmin,ymin,zmin,xmax,ymax,zmax, & 
     &                      lperiodicx,lperiodicy,lperiodicz, & 
     &                      coord(:,lb),bsize(:,lb),ndim, & 
     &                      lrefine(lb),lrefine_max,mort_neigh, & 
     &                      bnd_box(1,1,lb))


! Now get the possible neighbors of the current block^s parent
      if(parent(1,lb).gt.0) then
        pbsize(:) = bsize(:,lb)*2.               ! size of parent block
        ioff = mod(which_child(lb)-1,2)        ! coord for parent block
        joff = mod((which_child(lb)-1)/2,2)
        koff = mod((which_child(lb)-1)/4,2)
        if(ioff.eq.0) then
          pcoord(1) = bnd_box(2,1,lb)
        else
          pcoord(1) = bnd_box(1,1,lb)
        endif
        if(joff.eq.0) then
          pcoord(2) = bnd_box(2,2,lb)
        else
          pcoord(2) = bnd_box(1,2,lb)
        endif
        if(ndim.lt.2) pcoord(2) = coord(2,lb)
        if(koff.eq.0) then
          pcoord(3) = bnd_box(2,3,lb)
        else
          pcoord(3) = bnd_box(1,3,lb)
        endif
        if(ndim.lt.3) pcoord(3) = coord(3,lb)

        if(spherical_pm) then
! should try to fix this section so pbndbox gets values which are consistent
! with neighbors to the last digit
        pbndbox = bnd_box(:,:,lb)
        if(ioff.eq.0) then
          pbndbox(2,1) = pcoord(1) + pbsize(1)
        elseif(ioff.eq.1) then
          pbndbox(1,1) = pcoord(1) - pbsize(1)
        endif
        if(joff.eq.0) then
          pbndbox(2,2) = pcoord(2) + pbsize(2)
        elseif(joff.eq.1) then
          pbndbox(1,2) = pcoord(2) - pbsize(2)
        endif
        if(koff.eq.0) then
          pbndbox(2,3) = pcoord(3) + pbsize(3)
        elseif(koff.eq.1) then
          pbndbox(1,3) = pcoord(3) - pbsize(3)
        endif
        endif

        call morton_neighbors(xmin,ymin,zmin,xmax,ymax,zmax, & 
     &                        lperiodicx,lperiodicy,lperiodicz, & 
     &                        pcoord(:),pbsize(:),ndim, & 
     &                        lrefine(lb)-1,lrefine_max,pmort_neigh, & 
     &                        pbndbox(1,1))

      endif                      ! end of parent if test


!-------------

! If parent is a remote block then puts its address on the list of
! remote blocks which are required.
      if(parent(1,lb).gt.0.and.parent(2,lb).ne.mype) then
            istack = istack+1
#ifdef DEBUGZ
            if(istack.gt.npts_neigh) then
              write(*,*) 'morton_bnd : ', & 
     &                   'istack exceeds npts_neigh : ', & 
     &                   'possible solution - increase npts_neigh'
              call mpi_abort(amr_mpi_meshComm,ierrorcode,ierr)
            endif
#endif /* DEBUG */
            if(istack.gt.npts_neigh1) call expand_neigh_morts_mcomm
            neigh_morts(:,1,istack) = pmort_neigh(:,2,2,2)
            neigh_morts(6,2,istack) = lrefine(lb)-1
            neigh_morts(6,3,istack) = 14      ! marks as a full block request
#ifdef DEBUG
      write(*,*) 'parent pmort_neigh(6,2,2,2) ',pmort_neigh(6,2,2,2), & 
     &       ' istack ',istack,' pe ',mype,' parent ',parent(:,lb) & 
     &  ,' of block ',lb
#endif /* DEBUG */
      endif                     ! end of parent if test

!-------------

! Now start to build the array neigh_morts which is a list of possible
! remote blocks which will be needed.

! First add any neighbors of the current block, eliminating
! any local blocks on the list


      do k = 2-k3d,2+k3d
      do j = 2-k2d,2+k2d
      do i = 1,3
        if(i.ne.2.or.j.ne.2.or.k.ne.2) then

! if neighbor block exists at this refinement level
        if(mort_neigh(6,i,j,k).gt.-1) then

          lremote = is_remote(mort_neigh(1:6,i,j,k),lrefine(lb),mype)

          if(lremote) then 
            istack = istack+1
            if(istack.gt.npts_neigh1) call expand_neigh_morts_mcomm
            neigh_morts(:,1,istack) = mort_neigh(:,i,j,k)
            neigh_morts(6,2,istack) = lrefine(lb)
! compute message type - note this index is computed to reflect the part
! of the remote block to be acquired, not the part of the local blocks
! guardcells which will be filled.
            neigh_morts(6,3,istack) = (4-i)+((4-j)-1)*3+((4-k)-1)*9
            if(nguarda.gt.nmax_lays) neigh_morts(6,3,istack) = 14
#ifdef DEBUG 
            write(*,*) 'mpi_mort_comm pe ',mype, & 
     &             ' blk ',lb,' ijk ',i,j,k, & 
     &             ' neigh_morts ', & 
     &              neigh_morts(:,:,istack),' istack ',istack
#endif /* DEBUG */
          else
#ifdef DEBUG 
            write(*,*) 'pe ',mype,' blk ',lb,' ijk ',i,j,k,' local'
#endif /* DEBUG */
          endif                         ! if(lremote)
        endif                           ! if(mort_neigh(i,j,k).gt.-1)

        if (advance_all_levels) then


! Now consider required neighbors of the current blocks parent
        if(parent(1,lb).gt.0) then

! If interpolation order or number of guard cells is small enough, we can limit
! this list to neighbors of the parent relevant to the child block under consideration.
! Otherwise we consider all parent neighbors.
          ioff = mod(which_child(lb)-1,2) 
          joff = mod((which_child(lb)-1)/2,2)
          koff = mod((which_child(lb)-1)/4,2)
          ii = i
          jj = j
          kk = k
          if(interp_max_order.lt.nmax_lays) then
            if(ioff.eq.0.and.i.eq.3) ii = 2
            if(ioff.eq.1.and.i.eq.1) ii = 2
            if(joff.eq.0.and.j.eq.3) jj = 2
            if(joff.eq.1.and.j.eq.1) jj = 2
            if(koff.eq.0.and.k.eq.3) kk = 2
            if(koff.eq.1.and.k.eq.1) kk = 2
          endif
          lremote = .false.
          if (pmort_neigh(6,i,j,k) > -1) then

          if( morton_less_than(pmort_neigh(1:6,ii,jj,kk), & 
     &                         morton_limits(1:6,1,1,mype+1)) ) & 
     &         lremote = .true.
          if( morton_greater_than(pmort_neigh(1:6,ii,jj,kk), & 
     &                            morton_limits(1:6,1,2,mype+1)) ) & 
     &         lremote = .true.
          if( (morton_equal(pmort_neigh(1:6,ii,jj,kk), & 
     &                      morton_limits(1:6,1,1,mype+1))) & 
     &             .and. & 
     &        (lrefine(lb)-1.lt.morton_limits(6,2,1,mype+1)) ) & 
     &         lremote = .true.
          if( (morton_equal(pmort_neigh(1:6,ii,jj,kk), & 
     &                      morton_limits(1:6,1,2,mype+1))) & 
     &             .and. & 
     &        (lrefine(lb)-1.gt.morton_limits(6,2,2,mype+1)) ) & 
     &         lremote = .true.


          if(lremote) then 
            istack = istack+1
            if(istack.gt.npts_neigh1) call expand_neigh_morts_mcomm
            neigh_morts(:,1,istack) = pmort_neigh(:,ii,jj,kk)
            neigh_morts(6,2,istack) = lrefine(lb)-1
            neigh_morts(6,3,istack) = (4-ii)+((4-jj)-1)*3+((4-kk)-1)*9
            if(nguarda.gt.nmax_lays) neigh_morts(6,3,istack) = 14
#ifdef DEBUG 
            write(*,*) 'pe ',mype,' blk ',lb,' iijjkk ',ii,jj,kk, & 
     &             ' neigh_morts ', & 
     &              neigh_morts(:,:,istack),' istack ',istack
#endif /* DEBUG  */
          else
#ifdef DEBUG 
            write(*,*) 'pe ',mype,' blk ',lb,' iijjkk ', & 
     &         ii,jj,kk,' parent neigh local'
#endif /* DEBUG  */
          endif                      ! if(lremote)
          endif                      ! if (pmort_neigh > -1)
        endif                        ! if(parent(1,lb).gt.0)

      end if                    ! advance_all_levels

      endif                     ! if(i.ne.2.or.j.ne.2.or.k.ne.2)
      enddo
      enddo
      enddo

      enddo                     ! end loop over blocks


!--------------------------------------------------
      if(istack.gt.0) then
!--------------------------------------------------

      call compress_list(neigh_morts, & 
     &                   istack,no_of_remote_neighs,mype, & 
     &                   nprocs,l_on_pe)
      istack = no_of_remote_neighs

      end if

!--------------------------------------------------
       if(istack.gt.0) then
!--------------------------------------------------
!
! Step 5.
! Construct a list of all processors from which the local processor should
! request morton number information.


! non-zero elements of COMMATRIX define which processor pairs need to 
! exchange morton number lists.

        do i = 1,no_of_remote_neighs
          i_pe = 1
          j_pe = -1
          do while(  & 
     &       ( morton_greater_than(neigh_morts(1:6,1,i), & 
     &                             morton_limits(1:6,1,2,i_pe)) & 
     &                               .or. & 
     &         (morton_equal(neigh_morts(1:6,1,i), & 
     &                       morton_limits(1:6,1,2,i_pe)).and. & 
     &          neigh_morts(6,2,i).gt.morton_limits(6,2,2,i_pe)  )  ) & 
     &          .and. (i_pe.le.nprocs) & 
     &            )
             i_pe = i_pe + 1
             if (i_pe > nprocs) exit
          enddo
          if(i_pe.le.nprocs) j_pe = i_pe
!
! If block has been located then update commatrix
          if(j_pe.ne.-1)  & 
     &      commatrix_recv(j_pe) =  commatrix_recv(j_pe) + 1

        enddo

#ifdef DEBUG
        write(*,*) 'pe ',mype,' commatrix bef gather ', & 
     &             commatrix_recv(1:nprocs)
#endif /* DEBUG  */


! record the number of processors which will communicate with the
! local processor.
       no_of_comms_to_send = 0
       kstack = 0
       do i = 1,nprocs
         no_of_comms_to_send = no_of_comms_to_send + & 
     &                          min( 1, commatrix_recv(i) )
         if(commatrix_recv(i).gt.0) then
           kstack = kstack+1
           pe_source(kstack) = i
         endif
       enddo
#ifdef DEBUG
       write(*,*) 'pe ',mype,' no_of_comms_to_send ', & 
     &           no_of_comms_to_send
#endif /* DEBUG  */

!--------------------------------------------------
       endif                     ! end of istack if test
!--------------------------------------------------
!
! Step 6.
! provide the complete COMMATRIX to all processors

      call MPI_AlltoAll (commatrix_recv,       1,MPI_INTEGER, & 
     &                   commatrix_send,       1,MPI_INTEGER, & 
     &                   amr_mpi_meshComm,ierror)

#ifdef DEBUG
        write(*,*) 'pe ',mype,' commatrix ', & 
     &             commatrix_send(1:nprocs)

        write(*,'(" ")')
        write(*,'(" COMMUNICATION MATRIX1: mort_comm")')
        write(*,'(" ")')
        write(*,'("pe  ",i3," commatrix_send ", & 
     &       2i3)') mype,(commatrix_send(i),i=1,nprocs)
        write(*,'(" ")')
        Call MPI_BARRIER(amr_mpi_meshComm, ierr)

#endif /* DEBUG  */
!--------------------------------------------------
!
! Step 7.
! Compute the maximum amount of morton information which any processor
! is going to receive.


       max_no_to_be_received = 0

       iprocs = 0
       do j = 1,nprocs
          iprocs = iprocs + min(1,commatrix_recv(j))
       enddo
       max_no_to_be_received = max(1,iprocs)

#ifdef DEBUG
       write(*,*) 'pe ',mype,' max_no_to_be_received ', & 
     &           max_no_to_be_received
#endif /* DEBUG  */


!--------------------------------------------------
!
! Step 8.

       call MPI_ALLREDUCE(lnblocks,  & 
     &                    max_no_of_blocks, & 
     &                    1, & 
     &                    MPI_INTEGER, & 
     &                    MPI_MAX, & 
     &                    amr_mpi_meshComm, & 
     &                    ierror)

! Dynamically allocate memory to store the remote morton information.

       if(allocated(r_mortonbnd)) deallocate(r_mortonbnd)
       allocate( r_mortonbnd(6,3,max_no_of_blocks, & 
     &           max(1,max_no_to_be_received) ), & 
     &           stat = allocation_status)
       if(allocation_status > 0) then
          write(*,*) 'morton_bnd : allocation error'
          call mpi_abort(amr_mpi_meshComm,ierrorcode,ierr)
       endif


!--------------------------------------------------

       if(allocated(recvrequest)) deallocate( recvrequest )
       allocate ( recvrequest(nprocs) )

       if(allocated(recvstatus)) deallocate( recvstatus )
       allocate ( recvstatus(MPI_STATUS_SIZE,nprocs) )

!
! Step 9.
! Exchange morton information between processors.

      pe_source   = -1
      isize = 3*max_no_of_blocks*6
      k = 0
      r_mortonbnd = -1

      do i = 1,nprocs
         isrc = i-1
         idest= mype
#ifdef PM_UNIQUE_MPI_TAGS
         itag = isrc*nprocs + idest+1 + tag_offset
#else
         itag = tag_offset
#endif

                                ! receive to pe=j
         if((commatrix_recv(i).gt.0)) then
            k = k+1
            pe_source(k) = isrc+1
            call Mpi_Irecv(r_mortonbnd(1,1,1,k),isize,MPI_INTEGER, & 
     &           isrc ,itag,amr_mpi_meshComm,recvrequest(k),ierr)
         endif
      enddo

      ll = 0
      do j = 1,nprocs
          isrc = mype
          idest= j-1
#ifdef PM_UNIQUE_MPI_TAGS
          itag = isrc*nprocs + idest+1 + tag_offset
#else
          itag = tag_offset
#endif
                                 ! send from mype=i
          if(commatrix_send(j).gt.0) then
             ll = ll+1
             call MPI_Ssend(mortonbnd(1,1,1),isize,MPI_INTEGER, & 
     &            idest,itag,amr_mpi_meshComm,ierr)
          endif
      enddo

      no_of_mortonbnds_received = k

#ifdef PM_UNIQUE_MPI_TAGS
      tag_offset = (nprocs-1)*nprocs + nprocs + tag_offset
#endif

      if(k.gt.0) & 
     &    call MPI_Waitall(k,recvrequest,recvstatus, & 
     &                     ierror)


#ifdef DEBUG
      write(*,*) 'pe ',mype,' no_of_mortonbnds_received ', & 
     &          no_of_mortonbnds_received
      write(*,*) 'pe ',mype,' r_mortonbnd(6,:,1:15,1) ', & 
     &          r_mortonbnd(6,:,1:15,1)
      call amr_flush(6)
      Call MPI_BARRIER(amr_mpi_meshComm, ierr)
#endif /* DEBUG  */
        
!--------------------------------------------------
!
! Step 10.
! Loop over this processor^s list of required neighbor blocks,
! identifying their remote location from the morton information received
! in step 9.


        do i = 1,no_of_remote_neighs
          i_pe = 1
          j_pe = -1
          do while(  & 
     &      (  morton_greater_than(neigh_morts(1:6,1,i), & 
     &                             morton_limits(1:6,1,2,i_pe)) & 
     &                             .or. & 
     &        (morton_equal(neigh_morts(1:6,1,i), & 
     &                      morton_limits(1:6,1,2,i_pe)).and. & 
     &         neigh_morts(6,2,i).gt.morton_limits(6,2,2,i_pe)  )  ) & 
     &         .and. (i_pe.le.nprocs) & 
     &            )
            i_pe = i_pe + 1
            if (i_pe > nprocs) exit
          enddo
          if(i_pe.le.nprocs) j_pe = i_pe

          rem_block = -1
          rem_pe = j_pe

          kk = -1
          do k=1,no_of_mortonbnds_received
            if(pe_source(k).eq.rem_pe) kk = k 
          enddo
          if(kk.gt.0) then
          do j=1,max_no_of_blocks
            if( morton_equal(r_mortonbnd(1:6,1,j,kk), & 
     &                       neigh_morts(1:6,1,i)) .and. & 
     &          r_mortonbnd(6,2,j,kk).eq.neigh_morts(6,2,i) ) & 
     &          rem_block = j
          enddo
          endif
          if(rem_block.eq.-1) rem_pe = -1

#ifdef DEBUG 
          write(*,*) 'pe ',mype,' neigh i ',i,' rem_pe ', & 
     &            rem_pe,' kk ',kk,' rem_block ',rem_block
#endif /* DEBUG  */

! neigh_morts(1:2,no_of_remote_neighs) is now being used to store 
! the remote addresses of the required neighbors.
! Here proc nos. run from 1 to nprocs.

          neigh_morts(:,1,i) = rem_block
          neigh_morts(:,2,i) = rem_pe

#ifdef DEBUG 
          write(*,*) 'pe ',mype,' neigh i ',i,' address ', & 
     &            neigh_morts(:,:,i)
#endif /* DEBUG  */
        enddo

!--------------------------------------------------
!
! Step 11.
! Check for any non-existent blocks in the neigh_morts list
! and remove them. Then reset commatrix.


      if(allocated(indx)) deallocate(indx)
      allocate(indx(no_of_remote_neighs))


      indx = 0
      jstack = 0

      do i=1,no_of_remote_neighs
        if(neigh_morts(6,1,i).gt.-1) then
#ifdef DEBUG 
          write(*,*) 'pe ',mype,' stack entry ',neigh_morts(6,1,i), & 
     &     ' does exists - not to be removed '
#endif /* DEBUG  */
          jstack = jstack+1
          indx(jstack) = i
        endif
      enddo
      do j=1,jstack
        neigh_morts(6,:,j) = neigh_morts(6,:,indx(j))
#ifdef DEBUG 
        write(*,*) 'pe ',mype,' remaining stack entry ',j, & 
     & ' neigh_morts(:,j) ',neigh_morts(:,:,j)
#endif /* DEBUG  */
      enddo
      if(no_of_remote_neighs.gt.jstack) & 
     &      neigh_morts(6,:,jstack+1:no_of_remote_neighs) = -1
#ifdef DEBUG 
      write(*,*) 'pe ',mype,' removed stack items ',jstack+1, & 
     &       ' to ',no_of_remote_neighs
#endif /* DEBUG  */
      istack = jstack
      no_of_remote_neighs = istack

!--------------------------------------------------
! Step 12.
! Reconstruct commatrix.


! non-zero elements of COMMATRIX define which processor pairs need to 
! exchange morton number lists. 
        commatrix_send = 0
        commatrix_recv = 0
        do i = 1,no_of_remote_neighs
          i_pe = neigh_morts(6,2,i)
          commatrix_recv(i_pe) =  commatrix_recv(i_pe) + 1
        enddo

!
! Eliminate any r_mortonbnds layers which are no longer required.
        jstack = 0
        do i = 1,no_of_comms_to_send
          i_pe = pe_source(i)
          if(commatrix_recv(i_pe).gt.0) then
            jstack = jstack+1
            indx(jstack) = i
          endif
        enddo
        do j=1,jstack
          r_mortonbnd(:,:,:,j) = r_mortonbnd(:,:,:,indx(j))
        enddo
        no_of_mortonbnds_received = jstack            
#ifdef DEBUG
      write(*,*) 'pe ',mype,' revised no_of_mortonbnds_received ', & 
     &          no_of_mortonbnds_received
#endif /* DEBUG  */

! record the number of processors which will communicate with the
! local processor.
       pe_source = -1
       no_of_comms_to_send = 0
       kstack = 0
       do i = 1,nprocs
         no_of_comms_to_send = no_of_comms_to_send + & 
     &                          min( 1, commatrix_recv(i) )
         if(commatrix_recv(i).gt.0) then
           kstack = kstack+1
           pe_source(kstack) = i
         endif
       enddo
#ifdef DEBUG
       write(*,*) 'pe ',mype,' no_of_comms_to_send ', & 
     &           no_of_comms_to_send
#endif /* DEBUG  */

!--------------------------------------------------
!
! Step 13.
! Repeat Step 6.
! provide the complete COMMATRIX to all processors

      call MPI_AlltoAll (commatrix_recv,       1,MPI_INTEGER, & 
     &                   commatrix_send,       1,MPI_INTEGER, & 
     &                   amr_mpi_meshComm,ierror)

#ifdef DEBUG
        write(*,*) 'pe ',mype,' commatrix ', & 
     &             commatrix_send(1:nprocs)

        write(*,'(" ")')
        write(*,'(" COMMUNICATION MATRIX2: mort_comm")')
        write(*,'(" ")')
        write(*,'("pe  ",i3," commatrix_send ", & 
     &       2i3)') mype,(commatrix_send(i),i=1,nprocs)
        write(*,'(" ")')
        Call MPI_BARRIER(amr_mpi_meshComm, ierr)
#endif /* DEBUG  */


!--------------------------------------------------
!
! Step 20.
! Deallocate any memory which was dynamically allocated for local use in this
! routine.

       if(allocated(recvrequest)) deallocate( recvrequest )
       if(allocated(recvstatus)) deallocate( recvstatus )


!--------------------------------------------------
      deallocate(neigh_morts)
      deallocate(indx)

#ifdef DEBUG_FLOW_TRACE
      write(*,*) 'pe ',mype,' exiting mpi_mort_comm_for_surrblks'
#endif /* DEBUG_FLOW_TRACE */

      return

      contains
        subroutine expand_neigh_morts_mcomm

              if(allocated(tneigh_morts)) deallocate(tneigh_morts)
              allocate(tneigh_morts(6,3,npts_neigh2))
              tneigh_morts(:,:,:istack-1) = neigh_morts(:,:,:istack-1)
              npts_neigh1 = npts_neigh1 + 3000
              npts_neigh2 = npts_neigh2 + 3000
              deallocate(neigh_morts)
              allocate(neigh_morts(6,3,npts_neigh2))
              neigh_morts(:,:,:istack-1) = tneigh_morts(:,:,:istack-1)
              deallocate(tneigh_morts)

        end subroutine expand_neigh_morts_mcomm

      end subroutine mpi_mort_comm_for_surrblks



