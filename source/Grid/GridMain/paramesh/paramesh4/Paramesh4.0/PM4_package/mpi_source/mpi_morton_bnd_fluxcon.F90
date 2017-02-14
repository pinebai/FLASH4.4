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

      subroutine mpi_morton_bnd_fluxcon( & 
     &                        mype,nprocs,tag_offset)


!------------------------------------------------------------------------
!
! This routine calculates the morton number for each block on mype.
! It stores the result along with the refinement level of each block into
! the array mortonbnd, and distributes this array among all processors.
!
!
! Written :     Peter MacNeice  and Michael Gehmeyr          July 2000
!------------------------------------------------------------------------
!
! Arguments:
!      mype           rank of local processor
!
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use timings
      use mpi_morton
      use constants
      use paramesh_comm_data

      use paramesh_mpi_interfaces, only : mpi_amr_write_flux_comm, & 
     &                                    compress_list, & 
     &                                    morton_neighbors

      implicit none

      include 'mpif.h'

      integer, intent(in)    ::  mype,nprocs
      integer, intent(inout) ::  tag_offset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables

      real    :: eps,accuracy
      real    :: xmin, ymin, zmin
      real    :: xmax, ymax, zmax

      character*40,save :: c_routine = 'mpi_morton_bnd_flxcon'

      integer :: lb,i,j,ipe,j00
      integer :: level,jstack,iprocs
      integer :: morton(6)
      integer :: mort_neigh(6,3,3,3)
!     integer :: neigh_morts(6,3,npts_neigh),indx(npts_neigh)
      integer,dimension (:,:,:),allocatable:: neigh_morts
      integer,dimension (:,:,:),allocatable:: tneigh_morts
      integer,dimension (:)    ,allocatable:: indx
      integer :: istart,iend
      integer :: i_pe,j_pe,rem_block,rem_pe
      integer :: no_of_comm_procs
      integer :: ierrorcode,ierr,allocation_status,ierror
      integer :: l_nodetype(3,3,3)
      integer :: no_of_remote_neighs
      integer :: max_no_to_be_received
      integer :: no_of_comms_to_send
      integer :: max_no_of_blocks
      integer :: no_of_comms_to_receive
      integer :: istack, k, itemp, kstack, isize, isrc, idest
      integer :: itag, ll, kk, ij, ie, if1, jf1, kf1, if2, jf2, kf2
      integer :: nodetype_1, nodetype_2, jj, jp, ip, ii
      integer,dimension (:),  allocatable :: recvrequest
      integer,dimension (:,:),allocatable :: recvstatus

      logical :: llll
      logical :: lswap,lfound
      logical :: is_found 
      logical :: morton_greater_than
      logical,external :: morton_equal
      logical :: morton_less_than
      logical,save :: l_on_pe = .true.
!      logical,save :: l_on_pe = .false.


      double precision :: time1
      double precision :: time2
      double precision :: time3
      double precision :: time4
      integer          :: i_timer = 4        ! mpi_morton_bnd_fluxcon

      integer          :: npts_neigh1,npts_neigh2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef DEBUG_FLOW_TRACE
      write(*,*) 'entered mpi_morton_bnd_fluxcon'
#endif /* DEBUG_FLOW_TRACE */

      accuracy = 100./10.**precision(accuracy)
      eps = accuracy

!
! Step 1.

      if (timing_mpi) then
         time1 = mpi_wtime()
         time2 = mpi_wtime()
      endif

      npts_neigh1 = npts_neigh
      npts_neigh2 = npts_neigh+100
      allocate(neigh_morts(6,3,npts_neigh2))

!
!
! This routine assumes that the grid blocks are ordered by morton
! number and that any blocks with different refinement levels but
! the same morton number are ordered from coarse to fine.



!--------------------------------------------------

! Compute xmin,ymin,zmin,xmax,ymax,zmax or get them from storage
      xmin = grid_xmin
      ymin = grid_ymin
      zmin = grid_zmin
      xmax = grid_xmax
      ymax = grid_ymax
      zmax = grid_zmax

! Initializations
      no_of_comm_procs = 0
      no_of_remote_neighs = 0
      max_no_to_be_received = 0
      max_no_to_send = 0
      commatrix_send = 0
      commatrix_recv = 0
      pe_source = -1
      pe_destination = -1
      no_of_comms_to_send = 0

      if (timing_mpi) then
         time3 = mpi_wtime()
      endif


!FLASH edit 28th January 2008.  The valgrind memcheck tool complains that
!a branch based on the value of neigh_morts is flawed because 
!neigh_morts is not initialised.  As such, neigh_morts is now initialised 
!to -1.  Please note that this intialisation statement was previously here, 
!but was commented out, perhaps so for performance reasons.
!-------------------------------------------------------------------
      neigh_morts = -1
!-------------------------------------------------------------------


      if (timing_mpi) then
              timer_mpi_morton_bnd(i_timer,2) =   & 
     &                            timer_mpi_morton_bnd(i_timer,2) & 
     &                          + mpi_wtime() - time3
              timer_mpi_morton_bnd(i_timer,1) =   & 
     &                            timer_mpi_morton_bnd(i_timer,1) & 
     &                          + mpi_wtime() - time2
      endif

!--------------------------------------------------

      if (timing_mpi) then
         time2 = mpi_wtime()
         time3 = mpi_wtime()
         time4 = mpi_wtime()
      endif
!
! We want to construct a list of neighbors of leaf blocks which
! are refined. If we compute the morton number for a neighbor
! of a leaf block, it will be the same as the morton number of
! its leftmost child. Therefore if we look for the (mort,level)
! combination for this lower left child, finding it tells us that
! the leaf blocks neighbor is refined, and so we should add the
! neighbors (mort,level) to the list.

! Construct a list of the bottom left child of all neighbors of
! parents of leaf blocks. This list should include both remote
! and local children.

      istack = 0

#ifdef DEBUG
      call amr_flush(6)
      Call MPI_BARRIER(amr_mpi_meshComm, ierr)
      write(*,*) 'xmin,ymin,zmin,xmax,ymax,zmax ', & 
     & xmin,ymin,zmin,xmax,ymax,zmax
#endif /* DEBUG */

      do lb=1,lnblocks

      if(nodetype(lb).eq.1) then

!-------------

!
! Get morton numbers for neighbors of any leaf blocks
      mort_neigh  = -1
      call morton_neighbors(xmin,ymin,zmin,xmax,ymax,zmax, & 
     &                      lperiodicx,lperiodicy,lperiodicz, & 
     &                      coord(1,lb),bsize(1,lb),ndim, & 
     &                      lrefine(lb),lrefine_max,mort_neigh, & 
     &                      bnd_box(1,1,lb))
#ifdef DEBUG
      write(*,*) 'exited morton_neighbors for blk lb=',lb
#endif /* DEBUG */

!-------------

! Now cycle through this list of neighbors constructing the
! (morton no., refinement level) for their bottom left child.
      do k = 2-k3d,2+k3d
      do j = 2-k2d,2+k2d
      do i = 1,3
#ifdef DEBUG
      write(*,*) 'cycle ijk ',i,j,k
      write(*,*) 'cycle ijk mort_neigh ',mort_neigh(6,i,j,k)
#endif /* DEBUG */
        if(i.ne.2.or.j.ne.2.or.k.ne.2) then

          if(mort_neigh(6,i,j,k).gt.-1) then
            istack = istack+1
#ifdef DEBUGX
            if(istack.gt.npts_neigh) then
              write(*,*) c_routine,' : ', & 
     &                   'istack exceeds npts_neigh : ', & 
     &                   'possible solution - increase npts_neigh'
              call mpi_abort(amr_mpi_meshComm,ierrorcode,ierr)
            endif
#endif /* DEBUGX */
            if(istack.gt.npts_neigh1) then
              call expand_neigh_morts_flux
            endif
            neigh_morts(:,1,istack) = mort_neigh(:,i,j,k)
            neigh_morts(6,2,istack) = lrefine(lb)+1

            j00 = j

            if(spherical_pm) then
! if this block is a polar block then change the way j is applied
! in the formula for neigh_morts(6,3,istack)            
            if(lsingular_line) then
            if(abs(bnd_box(1,2,lb)).lt.eps.and.j.eq.1) then
             j00 = 3
            elseif(abs(bnd_box(2,2,lb)-pi).lt.eps.and.j.eq.3) then
             j00 = 1
            endif
            endif
            endif

            neigh_morts(6,3,istack) = (4-i)+((4-j00)-1)*3+((4-k)-1)*9
            if(nguard.gt.nmax_lays) neigh_morts(6,3,istack) = 14
#ifdef DEBUG
            write(*,*) 'flx:pe ',mype,' blk ',lb,' ijk ',i,j,k, & 
     &             ' neigh_morts ', & 
     &              neigh_morts(:,1,istack),' istack ',istack
#endif /* DEBUG */
          endif   

        endif
      enddo
      enddo
      enddo

!-------------

      if (timing_mpi) then
      timer_mpi_morton_bnd3(i_timer,4) =   & 
     &                            timer_mpi_morton_bnd3(i_timer,4) & 
     &                          + mpi_wtime() - time4
      time4 = mpi_wtime()
      endif

      endif

      enddo                     ! end loop over blocks


!--------------------------------------------------
      if(istack.gt.0) then
!--------------------------------------------------

! Calling compress_list here is a problem, since it will eliminate
! any on-proc neighbors, which we need. So we have added an additional
! argument which lets us keep these if we need them.
      call compress_list(neigh_morts,istack,no_of_remote_neighs,mype, & 
     &                   nprocs,l_on_pe)
      istack = no_of_remote_neighs

      end if

!--------------------------------------------------
       if(istack.gt.0) then
!--------------------------------------------------

!--------------------------------------------------
!
! Step 4.
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
!
! In this case some elements in the list might be local blocks. However
! we do not wish these to contribute to commatrix.
        commatrix_recv(mype+1) = 0

#ifdef DEBUGY
        write(*,*) 'pe ',mype,' commatrix bef gather '
!     .             ,commatrix(1:nprocs,1:nprocs)
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
! Step 5.
! provide the complete COMMATRIX to all processors

      call MPI_AlltoAll (commatrix_recv,       1,MPI_INTEGER, & 
     &                   commatrix_send,       1,MPI_INTEGER, & 
     &                   amr_mpi_meshComm,ierror)
       
#ifdef DEBUGY
        write(*,*) 'pe ',mype,' commatrix ', & 
     &             commatrix_send(1:nprocs)
 
         write(*,'(" ")')
         write(*,'(" COMMUNICATION MATRIX: m_bnd")')
         write(*,'(" ")')
         write(*,'("pe  ",i3," commatrix_send ", & 
     &       2i3)') mype,(commatrix_send(i),i=1,nprocs)
         write(*,'(" ")')
        Call MPI_BARRIER(amr_mpi_meshComm, ierr)

#endif /* DEBUG  */
!--------------------------------------------------
!
! Step 6.
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
! Step 7.
! Dynamically allocate memory to store the remote morton information.

       call MPI_ALLREDUCE(lnblocks, & 
     &                    max_no_of_blocks, & 
     &                    1, & 
     &                    MPI_INTEGER, & 
     &                    MPI_MAX, & 
     &                    amr_mpi_meshComm, & 
     &                    ierror)

       if(allocated(r_mortonbnd)) deallocate(r_mortonbnd)
       allocate( r_mortonbnd(6,3,max_no_of_blocks, & 
     &                       max(1,max_no_to_be_received) ), & 
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

!--------------------------------------------------
!
! Step 8.
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
          if(commatrix_recv(i).gt.0) then
             k = k+1
             pe_source(k) = isrc+1
             call Mpi_Irecv(r_mortonbnd(1,1,1,k),isize,MPI_INTEGER, & 
     &            isrc ,itag,amr_mpi_meshComm,recvrequest(k),ierr)
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
     &           idest,itag,amr_mpi_meshComm,ierr)
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
!      write(*,*) 'pe ',mype,' r_mortonbnd(:,1:15,1) ',
!     .          r_mortonbnd(:,1:15,1)
#endif /* DEBUG  */
        
!--------------------------------------------------
!
! Step 9.
! Loop over this processor^s list of required neighbor blocks,
! identifying whether they exist from the morton information received
! in step 8.

        do i = 1,no_of_remote_neighs
#ifdef DEBUG
      write(*,*) 'pe ',mype,' i=',i,' line A no_of_remote_neighs=', & 
     &           no_of_remote_neighs
#endif /* DEBUG  */

! first test whether this block can be found locally
!          write(*,*) 'i= ',i,' neigh_morts(:,1,i) ',neigh_morts(:,1,i)
!          write(*,*) 'i= ',i,' neigh_morts(6,2,i) ',neigh_morts(6,2,i)
!          write(*,*) 'i= ',i,' mortonbnd(:,1,510) ',mortonbnd(:,1,510)
!          write(*,*) 'i= ',i,' mortonbnd(6,2,510) ',mortonbnd(6,2,510)
!          write(*,*) 'i= ',i,' max_no_of_blocks ',max_no_of_blocks
          rem_block = -1
          rem_pe = mype+1
          do j=1,max_no_of_blocks
            if( morton_equal(mortonbnd(1:6,1,j),neigh_morts(1:6,1,i))  & 
     &                    .and. & 
     &          mortonbnd(6,2,j).eq.neigh_morts(6,2,i) ) & 
     &          rem_block = j
          enddo
! If not found locally, test to see if it is in any of the remote lists
! which have been received.
#ifdef DEBUG
      write(*,*) 'pe ',mype,' i=',i,' line B'
#endif /* DEBUG  */

          if(rem_block.eq.-1) then

            i_pe = 1
            j_pe = -1
            do while(  & 
     &        (  morton_greater_than(neigh_morts(1:6,1,i), & 
     &                               morton_limits(1:6,1,2,i_pe)) & 
     &                             .or. & 
     &          (morton_equal(neigh_morts(1:6,1,i), & 
     &                        morton_limits(1:6,1,2,i_pe)).and. & 
     &           neigh_morts(6,2,i).gt.morton_limits(6,2,2,i_pe)  )  ) & 
     &           .and. (i_pe.le.nprocs) & 
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
     &                         neigh_morts(1:6,1,i)) .and. & 
     &            r_mortonbnd(6,2,j,kk).eq.neigh_morts(6,2,i) ) then
                 rem_block = j

              endif
            enddo
            endif
            if(rem_block.eq.-1) rem_pe = -1

          else

          endif

#ifdef DEBUG
      write(*,*) 'pe ',mype,' i=',i,' line C'
#endif /* DEBUG  */
!
! If the child block exists then record the (mort no., level) pair 
! of it^s parent.

          if(rem_block.lt.0) then
            neigh_morts(:,:,i) = -1
          else
            neigh_morts(6,2,i) = neigh_morts(6,2,i)-1
          endif

! At this point neigh_morts has been modified to
! contain a list of (morton no, refinement level) pairs associated
! with the parents of any children which were successfully located
! during the first sweep.

        enddo

#ifdef DEBUG
      write(*,*) 'pe ',mype,' mortonfluxcon step 9 complete'
#endif /* DEBUG  */
        
!--------------------------------------------------
!
! Step 10.
! Check for any non-existent blocks in the neigh_morts list
! and remove them, compressing the list.

      if(allocated(indx)) deallocate(indx)
      allocate(indx(no_of_remote_neighs))


      indx = 0
      jstack = 0

      do i=1,no_of_remote_neighs
        if(neigh_morts(6,1,i).gt.-1) then
          jstack = jstack+1
          indx(jstack) = i
        endif
      enddo
      do j=1,jstack
        neigh_morts(:,:,j) = neigh_morts(:,:,indx(j))
      enddo
      neigh_morts(:,:,jstack+1:istack) = -1
      istack = jstack
      no_of_remote_neighs = istack

#ifdef DEBUG
      write(*,*) 'pe ',mype,' mortonfluxcon step 10 complete'
#endif /* DEBUG  */
!--------------------------------------------------

      if (nedgevar1 > 0) then

!--------------------------------------------------
!
! Step 11.
! If edge averaging is required then we need to identify leaf blocks
! which have refinement off their diagonal edges but not in the neighbors
! bounding that diagonal element. In this case we will have to make
! sure that the edge integral on the local block match those of the
! refined diagonal neighbor, if we wish to preserve a div B constraint.
!

      jstack = 0
      do lb = 1,lnblocks
        if(nodetype(lb).eq.1) then

! Compute the morton numbers of the children on diagonal edges

! Get morton numbers for neighbors of any leaf blocks
           mort_neigh  = -1
           call morton_neighbors(xmin,ymin,zmin,xmax,ymax,zmax, & 
     &                      lperiodicx,lperiodicy,lperiodicz, & 
     &                      coord(1,lb),bsize(1,lb),ndim, & 
     &                      lrefine(lb),lrefine_max,mort_neigh, & 
     &                      bnd_box(1,1,lb))

           l_nodetype = 1
           do k = 2-k3d,2+k3d
           do j = 2-k2d,2+k2d
           do i = 1,3
             if(i.ne.2.or.j.ne.2.or.k.ne.2) then
               lfound = .false.
               if(mort_neigh(6,i,j,k).gt.-1) then
                 morton(:) = mort_neigh(:,i,j,k)
                 level= lrefine(lb)
                 do ij = 1,no_of_remote_neighs
                   if( morton_equal(morton(1:6),neigh_morts(1:6,1,ij)) & 
     &                      .and. & 
     &                 level.eq.neigh_morts(6,2,ij) ) then
                     lfound = .true.
                   endif
                 enddo
               endif
               if(lfound) l_nodetype(i,j,k) = 2
             endif
           enddo
           enddo
           enddo

! Now cycle through this list of neighbors constructing the
! (morton no., refinement level) for their bottom left child.
      do ie = 1,12
        i = 2
        j = 2
        k = 2
        if1 = 2
        jf1 = 2
        kf1 = 2
        if2 = 2
        jf2 = 2
        kf2 = 2
        if(ie.eq.1) then
          i = 1
          j = 1
          if1 = 1
          jf2 = 1
        elseif(ie.eq.2) then
          i = 1
          j = 3
          if1 = 1
          jf2 = 3
        elseif(ie.eq.3) then
          i = 3
          j = 1
          if1 = 3
          jf2 = 1
        elseif(ie.eq.4) then
          i = 3
          j = 3
          if1 = 3
          jf2 = 3
        elseif(ie.eq.5) then
          j = 1
          k = 1
          jf1 = 1
          kf2 = 1
        elseif(ie.eq.6) then
          j = 3
          k = 1
          jf1 = 3
          kf2 = 1
        elseif(ie.eq.7) then
          j = 1
          k = 3
          jf1 = 1
          kf2 = 3
        elseif(ie.eq.8) then
          j = 3
          k = 3
          jf1 = 3
          kf2 = 3
        elseif(ie.eq.9) then
          i = 1
          k = 1
          if1 = 1
          kf2 = 1
        elseif(ie.eq.10) then
          i = 1
          k = 3
          if1 = 1
          kf2 = 3
        elseif(ie.eq.11) then
          i = 3
          k = 1
          if1 = 3
          kf2 = 1
        elseif(ie.eq.12) then
          i = 3
          k = 3
          if1 = 3
          kf2 = 3
        endif

           lfound = .false.
           if(mort_neigh(6,i,j,k).gt.-1) then
             morton(:) = mort_neigh(:,i,j,k)
             level= lrefine(lb)
! find (mort,level) in the list of diagonal edges
             do istack = 1,no_of_remote_neighs
               if( morton_equal(morton(1:6),neigh_morts(1:6,1,istack)) & 
     &                  .and. & 
     &             level.eq.neigh_morts(6,2,istack) ) then
                 lfound = .true.
               endif
             enddo
           endif
           if(lfound) then
! Are the neighbors which bound this edge refined?
             nodetype_1 = l_nodetype(if1,jf1,kf1)
             nodetype_2 = l_nodetype(if2,jf2,kf2)
             if(nodetype_1.eq.1.and.nodetype_2.eq.1) then
               jstack = jstack+1
               edge_mark(6,1,jstack) = ie
               edge_mark(6,2,jstack) = lb
               edge_mark(:,3,jstack) = morton(:)
               edge_mark(6,4,jstack) = level
             endif
           else
           endif

      enddo


        endif
      enddo
      no_of_diagonal_edges = jstack

      end if

#ifdef DEBUG
      write(*,*) 'pe ',mype,' mortonfluxcon step 11 complete'
#endif /* DEBUG  */

!--------------------------------------------------
!
! Step 12.
! Now that we have identified the edges which will need diagonal correction,
! we no longer need any local entries in the neigh_morts list.

!
! Mark any entries in the neigh_morts list which are on this processor
! for removal from the list.

          do i=1,no_of_remote_neighs

            if( & 
     &      (  morton_greater_than(neigh_morts(1:6,1,i), & 
     &                             morton_limits(1:6,1,1,mype+1)) & 
     &                             .or. & 
     &        (morton_equal(neigh_morts(1:6,1,i), & 
     &                      morton_limits(1:6,1,1,mype+1)).and. & 
     &         neigh_morts(6,2,i).ge.morton_limits(6,2,1,mype+1)  )  ) & 
     &    .and. & 
     &      (  morton_less_than(neigh_morts(1:6,1,i), & 
     &                          morton_limits(1:6,1,2,mype+1)) & 
     &                             .or. & 
     &        (morton_equal(neigh_morts(1:6,1,i), & 
     &                      morton_limits(1:6,1,2,mype+1)).and. & 
     &         neigh_morts(6,2,i).le.morton_limits(6,2,2,mype+1)  )  ) & 
     &      ) then
              neigh_morts(:,:,i) = -1
            endif

          enddo

#ifdef DEBUG
      write(*,*) 'pe ',mype,' mortonfluxcon step 12 complete'
#endif /* DEBUG  */
!--------------------------------------------------
!
! Step 13.
! Check for any non-existent blocks in the neigh_morts list
! and remove them.

      indx = 0
      jstack = 0
      do i=1,no_of_remote_neighs
        if(neigh_morts(6,1,i).gt.-1) then
          jstack = jstack+1
          indx(jstack) = i
        endif
      enddo
      do j=1,jstack
        neigh_morts(:,:,j) = neigh_morts(:,:,indx(j))
      enddo
      if(no_of_remote_neighs.gt.jstack) & 
     &      neigh_morts(:,:,jstack+1:istack) = -1
      istack = jstack
      no_of_remote_neighs = istack

!--------------------------------------------------
!
! If no such blocks were detected exit this routine.

!     if(no_of_remote_neighs.eq.0) return


!--------------------------------------------------
!
! The next stage is to locate the neighbors of the leaf
! blocks which were identified in the previous section.
! All these blocks will exist but they may be on or off
! processor.

!
! ReInitializations
      no_of_comm_procs = 0
      max_no_to_be_received = 0
      max_no_to_send = 0
      commatrix_recv = 0
      commatrix_send = 0
      pe_source = -1
      pe_destination = -1
!
#ifdef DEBUG
      write(*,*) 'pe ',mype,' mortonfluxcon step 13 complete'
#endif /* DEBUG  */
! Step 14.
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
         write(*,'(" COMMUNICATION MATRIX1: m_bnd")')
         write(*,'(" ")')
         write(*,'("pe  ",i3," commatrix_send ", & 
     &       2i3)') mype,(commatrix_send(i),i=1,nprocs)
         write(*,'(" ")')
        Call MPI_BARRIER(amr_mpi_meshComm, ierr)
#endif /* DEBUG  */
      if (timing_mpi) then
              timer_mpi_morton_bnd(i_timer,6) =   & 
     &                            timer_mpi_morton_bnd(i_timer,6) & 
     &                          + mpi_wtime() - time2
      endif

!--------------------------------------------------

      if (timing_mpi) then
         time2 = mpi_wtime()
      endif
!
! Step 7.
! Compute the maximum amount of morton information which any processor
! is going to receive.

       iprocs = 0
       do j = 1,nprocs
          iprocs = iprocs + min(1,commatrix_recv(j))
       enddo
       max_no_to_be_received = max(1,iprocs)

#ifdef DEBUG
       write(*,*) 'pe ',mype,' max_no_to_be_received ', & 
     &           max_no_to_be_received
#endif /* DEBUG  */

       if (timing_mpi) then
              timer_mpi_morton_bnd(i_timer,7) =   & 
     &                            timer_mpi_morton_bnd(i_timer,7) & 
     &                          + mpi_wtime() - time2
       endif

!--------------------------------------------------

      if (timing_mpi) then
         time2 = mpi_wtime()
      endif

!
! Step 14.
! Dynamically allocate memory to store the remote morton information.

       if(allocated(r_mortonbnd)) deallocate(r_mortonbnd)
       allocate( r_mortonbnd(6,3,max_no_of_blocks, & 
     &           max(1,max_no_to_be_received) ), & 
     &           stat = allocation_status)
       if(allocation_status > 0) then
          write(*,*) c_routine,' : allocation error'
          call mpi_abort(amr_mpi_meshComm,ierrorcode,ierr)
       endif


!--------------------------------------------------
!
! Step 15.
! Exchange morton information between processors.

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
#endif /* DEBUG  */
        

      if (timing_mpi) then
              timer_mpi_morton_bnd(i_timer,9) =   & 
     &                            timer_mpi_morton_bnd(i_timer,9) & 
     &                          + mpi_wtime() - time2
      endif

!--------------------------------------------------

      if (timing_mpi) then
         time2 = mpi_wtime()
      endif

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

            if(rem_block.lt.0.or.rem_pe.eq.mype+1) then
              neigh_morts(:,:,i) = -1
            else
              neigh_morts(:,1,i) = rem_block
              neigh_morts(:,2,i) = rem_pe
            endif

! At this point neigh_morts has been modified to
! contain the addresses of neighbors of local leaf
! blocks, provided the neighbors are parents and are located
! off processor.

        enddo


        if (nedgevar1 > 0) then
!
! Step 16a.
! Repeat step 16 but this time replace edge_mark(3:4,..) with
! the appropriate block address.

        do i = 1,no_of_diagonal_edges

          i_pe = 1
          j_pe = -1

          do while(  & 
     &      (  morton_greater_than(edge_mark(1:6,3,i), & 
     &                             morton_limits(1:6,1,2,i_pe)) & 
     &                             .or. & 
     &        (morton_equal(edge_mark(1:6,3,i), & 
     &                      morton_limits(1:6,1,2,i_pe)).and. & 
     &         edge_mark(6,4,i).gt.morton_limits(6,2,2,i_pe)  )  ) & 
     &         .and. (i_pe.le.nprocs) & 
     &            )
            i_pe = i_pe + 1
            if (i_pe > nprocs) exit
          enddo
          if(i_pe.le.nprocs) j_pe = i_pe

          rem_block = -1
          rem_pe = j_pe

          if(rem_pe.eq.mype+1) then

            do j=1,max_no_of_blocks
              if( morton_equal(mortonbnd(1:6,1,j), & 
     &                         edge_mark(1:6,3,i)) .and. & 
     &            mortonbnd(6,2,j).eq.edge_mark(6,4,i) ) & 
     &          rem_block = j
            enddo

          else

            kk = -1
            do k=1,no_of_mortonbnds_received
              if(pe_source(k).eq.rem_pe) kk = k 
            enddo
            if(kk.gt.0) then
            do j=1,max_no_of_blocks
              if( morton_equal(r_mortonbnd(1:6,1,j,kk), & 
     &                         edge_mark(1:6,3,i)) .and. & 
     &            r_mortonbnd(6,2,j,kk).eq.edge_mark(6,4,i) ) & 
     &            rem_block = j
            enddo
            endif

          endif
          if(rem_block.eq.-1) rem_pe = -1

#ifdef DEBUG 
          write(*,*) 'pe ',mype,' neigh i ',i,' rem_pe ', & 
     &            rem_pe,' kk ',kk,' rem_block ',rem_block
#endif /* DEBUG  */

            if(rem_block.lt.0) then
              edge_mark(:,3:4,i) = -1
            else
              edge_mark(6,3,i) = rem_block
              edge_mark(6,4,i) = rem_pe-1
            endif

        enddo
        end if

        if (timing_mpi) then
              timer_mpi_morton_bnd(i_timer,10) =   & 
     &                            timer_mpi_morton_bnd(i_timer,10) & 
     &                          + mpi_wtime() - time2
        endif

!--------------------------------------------------

      if (timing_mpi) then
         time2 = mpi_wtime()
      endif
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
     & ' neigh_morts(:,j) ',neigh_morts(6,:,j)
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

      if (timing_mpi) then
              timer_mpi_morton_bnd(i_timer,11) =   & 
     &                           timer_mpi_morton_bnd(i_timer,11) & 
     &                          + mpi_wtime() - time2
      endif

!--------------------------------------------------

      if (timing_mpi) then
         time2 = mpi_wtime()
      endif

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

       if (timing_mpi) then
              timer_mpi_morton_bnd(i_timer,12) =   & 
     &                            timer_mpi_morton_bnd(i_timer,12) & 
     &                          + mpi_wtime() - time2
       endif

!--------------------------------------------------

      if (timing_mpi) then
         time2 = mpi_wtime()
      endif

!
! Step 13.
! Repeat Step 6.
! provide the complete COMMATRIX to all processors

      call MPI_AlltoAll (commatrix_recv,       1,MPI_INTEGER, & 
     &                   commatrix_send,       1,MPI_INTEGER, & 
     &                   amr_mpi_meshComm,ierror)

#ifdef DEBUG
        write(*,*) 'pe ',mype,' commatrix ', & 
     &             commatrix_recv(1:nprocs)


         write(*,'(" ")')
         write(*,'(" COMMUNICATION MATRIX2: m_bnd")')
         write(*,'(" ")')
         write(*,'("pe  ",i3," commatrix_send ", & 
     &       2i3)') mype,(commatrix_send(i),i=1,nprocs)
         write(*,'(" ")')
        Call MPI_BARRIER(amr_mpi_meshComm, ierr)
#endif /* DEBUG  */

      if (timing_mpi) then
              timer_mpi_morton_bnd(i_timer,13) =   & 
     &                            timer_mpi_morton_bnd(i_timer,13) & 
     &                          + mpi_wtime() - time2
      endif

!--------------------------------------------------
      if (timing_mpi) then
         time2 = mpi_wtime()
      endif

! Step 14.
! record the number of processors to which the local processor
! will send messages.

       no_of_comms_to_receive = 0
       kstack = 0
       do i = 1,nprocs
         no_of_comms_to_receive = no_of_comms_to_receive + & 
     &                          min( 1, commatrix_send(i) )
         if(commatrix_send(i).gt.0) then
           kstack = kstack+1
           pe_destination(kstack) = i
         endif
       enddo
#ifdef DEBUG
       write(*,*) 'pe ',mype,' no_of_comms_to_receive ', & 
     &           no_of_comms_to_receive
#endif /* DEBUG  */

       if (timing_mpi) then
              timer_mpi_morton_bnd(i_timer,14) =   & 
     &                            timer_mpi_morton_bnd(i_timer,14) & 
     &                          + mpi_wtime() - time2
       end if

!--------------------------------------------------
       if (timing_mpi) then
          time2 = mpi_wtime()
       endif
!
! Step 15.
! Compute the maximum amount of morton information which any processor
! is going to receive.

       iprocs = 0
       do j = 1,nprocs
          iprocs = iprocs + min(1,commatrix_recv(j))
       enddo
       max_no_to_be_received = max(1,iprocs)

#ifdef DEBUGX
       write(*,*) 'pe ',mype,' max_no_to_be_received ', & 
     &           max_no_to_be_received
#endif /* DEBUGX  */

       if (timing_mpi) then
              timer_mpi_morton_bnd(i_timer,15) =   & 
     &                            timer_mpi_morton_bnd(i_timer,15) & 
     &                          + mpi_wtime() - time2
       endif
!--------------------------------------------------
       if (timing_mpi) then
          time2 = mpi_wtime()
       endif
!
! Step 16.
! Compute the maximum amount of information which any processor
! is going to receive.

       iprocs = 0
       do j = 1,nprocs
          iprocs = iprocs + min(1,commatrix_send(j))
       enddo
       max_no_to_send = max(1,iprocs)

#ifdef DEBUG
       write(*,*) 'pe ',mype,' max_no_to_send ', & 
     &           max_no_to_send
#endif /* DEBUG  */

       if (timing_mpi) then
              timer_mpi_morton_bnd(i_timer,16) =   & 
     &                            timer_mpi_morton_bnd(i_timer,16) & 
     &                          + mpi_wtime() - time2
       endif

!--------------------------------------------------
      if (timing_mpi) then
         time2 = mpi_wtime()
      endif

!
! Step 17.
! evaluate smallest guard block starting index over all pe
! store this into variable strt_buffer.

      last_buffer = maxblocks_alloc

      k = last_buffer
      do i=0,nprocs-1
      k = k - commatrix_recv(i+1)
      enddo
      strt_buffer = k + 1


      if (strt_buffer.le.lnblocks) then
        write(*,*)  & 
     &  'ERROR in ',c_routine,' : guard block starting index', & 
     &  strt_buffer,' not larger than lnblocks',lnblocks, & 
     &  ' processor no. ',mype,' maxblocks_alloc ', & 
     &  maxblocks_alloc
        call mpi_abort(amr_mpi_meshComm,ierrorcode,ierr)
      endif

      if (timing_mpi) then
              timer_mpi_morton_bnd(i_timer,17) =   & 
     &                            timer_mpi_morton_bnd(i_timer,17) & 
     &                          + mpi_wtime() - time2
      endif

!--------------------------------------------------
      if (timing_mpi) then
         time2 = mpi_wtime()
      endif
!
! Step 18.
! Dynamically allocate memory to store the lists of blocks to be
! sent and received.

      iprocs = max(maxval(commatrix_send),maxval(commatrix_recv))
      call MPI_ALLREDUCE(iprocs,  & 
     &                   largest_no_of_blocks, & 
     &                   1, & 
     &                   MPI_INTEGER, & 
     &                   MPI_MAX, & 
     &                   amr_mpi_meshComm, & 
     &                   ierror)
      
       if(allocated(to_be_sent)) deallocate(to_be_sent)
       if(allocated(to_be_received)) deallocate(to_be_received)
       allocate( to_be_sent(3, & 
     &                            max(1,largest_no_of_blocks), & 
     &                            max(1,max_no_to_send) ) )
       allocate( to_be_received(3, & 
     &                          max(1,largest_no_of_blocks), & 
     &                          max(1,max_no_to_be_received) ) )

       if (timing_mpi) then
              timer_mpi_morton_bnd(i_timer,18) =   & 
     &                            timer_mpi_morton_bnd(i_timer,18) & 
     &                          + mpi_wtime() - time2
       endif

!--------------------------------------------------
      if (timing_mpi) then
         time2 = mpi_wtime()
      endif
!
! Step 19.
! Construct arrays to_be_sent and to_be_received which contain
! the lists of blocks to be packaged.


        to_be_sent = -1
        to_be_received = -1
        laddress = 0

! First set up the array to_be_received on each processor
        if(no_of_remote_neighs.gt.0) then

          jj = 0
          do jp = 1,no_of_mortonbnds_received
            ip = pe_source(jp)
            if(commatrix_recv(ip).gt.0) then   ! this is a needless check
              do ii = 1,commatrix_recv(ip)
                jj = jj+1
                if(neigh_morts(6,2,jj).eq.ip) then
                  if(ii.gt.largest_no_of_blocks) then
          write(*,*) 'pe ',mype,' ii too large ',ii
                  endif
                  to_be_received(:,ii,jp) = neigh_morts(6,:,jj)
                endif
              enddo 
            endif
          enddo

          laddress(1,strt_buffer:strt_buffer+jj-1) = & 
     &          neigh_morts(6,1,1:jj)
          laddress(2,strt_buffer:strt_buffer+jj-1) = & 
     &          neigh_morts(6,2,1:jj)-1

        endif

#ifdef DEBUGX
        do jp = 1,no_of_mortonbnds_received
        write(*,*) 'pe ',mype,' jreceive ',jp,' to_be_received ', & 
     &    to_be_received(:,:,jp)
        enddo
#endif /* DEBUGX  */


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now exchange info in to_be_received with the sending processors
! to construct the equivalent to_be_sent arrays

! Post receives 
        isize = 3*largest_no_of_blocks
        kk = 0
        do i = 1,nprocs

          isrc = i-1
          idest= mype
#ifdef PM_UNIQUE_MPI_TAGS
          itag = isrc*nprocs + idest+1 + tag_offset
#else
          itag = tag_offset
#endif

                                 ! receive to pe=j
          if(commatrix_send(i).gt.0) then
            kk = kk+1
            call Mpi_Irecv(to_be_sent(1,1,kk),isize, & 
     &                     MPI_INTEGER,isrc ,itag,amr_mpi_meshComm, & 
     &                     recvrequest(kk),ierr)
          endif
        enddo

! Post sends

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
          if(commatrix_recv(j).gt.0) then
            ll = ll+1
            call MPI_Ssend(to_be_received(1,1,ll),isize,MPI_INTEGER, & 
     &           idest,itag,amr_mpi_meshComm,ierr)
          endif
        enddo

#ifdef PM_UNIQUE_MPI_TAGS
        tag_offset = (nprocs-1)*nprocs + nprocs + tag_offset
#endif

        if(kk.gt.0) & 
     &    call MPI_Waitall(kk,recvrequest,recvstatus, & 
     &                     ierror)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef DEBUG
        do jp = 1,kk
        write(*,*) 'pe ',mype,' jsend ',jp,' to_be_sent ', & 
     &    to_be_sent(:,:,jp)
        enddo
#endif /* DEBUG */

        if (timing_mpi) then
              timer_mpi_morton_bnd(i_timer,19) =   & 
     &                            timer_mpi_morton_bnd(i_timer,19) & 
     &                          + mpi_wtime() - time2
        endif

!--------------------------------------------------
      if (timing_mpi) then
         time2 = mpi_wtime()
      endif
!
! Step 20.
! Deallocate any memory which was dynamically allocated for local use in this
! routine.

       if(allocated(recvrequest)) deallocate( recvrequest )
       if(allocated(recvstatus)) deallocate( recvstatus )


!--------------------------------------------------

! Mark morton data up to date
       morton_limits_set = .true.


! Store communication info for future use
       call mpi_amr_write_flux_comm(nprocs)



#ifdef DEBUG
      write(*,*) 'pe ',mype,' exiting ',c_routine
#endif /* DEBUG */

      if (timing_mpi) then
              timer_mpi_morton_bnd(i_timer,20) =   & 
     &                            timer_mpi_morton_bnd(i_timer,20) & 
     &                          + mpi_wtime() - time2
              timer_mpi_morton_bnd(i_timer,0) =   & 
     &                            timer_mpi_morton_bnd(i_timer,0) & 
     &                          + mpi_wtime() - time1
      endif


      deallocate(neigh_morts)
      deallocate(indx)

#ifdef DEBUG_FLOW_TRACE
      write(*,*) 'exiting mpi_morton_bnd_fluxcon'
#endif /* DEBUG_FLOW_TRACE */

      return

      contains
        subroutine expand_neigh_morts_flux

              if(allocated(tneigh_morts)) deallocate(tneigh_morts)
              allocate(tneigh_morts(6,3,npts_neigh2))
              tneigh_morts(:,:,:istack-1) = neigh_morts(:,:,:istack-1)
              npts_neigh1 = npts_neigh1 + 3000
              npts_neigh2 = npts_neigh2 + 3000
              deallocate(neigh_morts)
              allocate(neigh_morts(6,3,npts_neigh2))
              neigh_morts(:,:,:istack-1) = tneigh_morts(:,:,:istack-1)
              deallocate(tneigh_morts)

        end subroutine expand_neigh_morts_flux

      end subroutine mpi_morton_bnd_fluxcon
