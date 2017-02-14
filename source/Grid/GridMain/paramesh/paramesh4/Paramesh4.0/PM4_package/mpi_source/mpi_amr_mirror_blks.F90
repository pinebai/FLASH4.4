!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/amr_mirror_blks
!! NAME
!!
!!   amr_mirror_blks
!!
!! SYNOPSIS
!!
!!   call amr_mirror_blks(mype, nprocs. mirror_blks, ix, iy, iz)
!!
!!   call amr_mirror_blks(integer, integer, integer, integer, integer, integer)
!!
!! ARGUMENTS
!!   
!!   integer, intent(in) :: mype  
!!     The calling processor.
!!
!!   integer, intent(in) :: nprocs
!!     The number for processors used.
!!
!!   integer, intent(inout) :: mirror_blks
!!     The list of mirror blocks found.
!!
!!   integer, intent(in) :: ix, iy, iz
!!     ???
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!   mpif.h
!!
!! USES
!! 
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   mpi_morton
!!   paramesh_interfaces   
!!   paramesh_mpi_interfaces
!!
!! CALLS
!! 
!!   morton_sort
!!   morton_number
!!    
!! RETURNS
!!
!!  The list of mirror blocks is returned in the array 'mirror_blks'.
!!
!! DESCRIPTION
!!
!!  This routine finds 'mirror' blocks (if they exist).  A mirror block which is
!!  a block symmetrically placed relative to some slice through the computational
!!  domain.  If the grid is not symmetric, an error message is issued.
!!  This routine is useful for debugging when the problem being run is
!!  supposed to maintain symmetry.
!!
!! AUTHORS
!!
!!   Peter MacNeice (December 2002).
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

!#define DEBUG
!#define DEBUGX
!#define DEBUGY

      subroutine amr_mirror_blks(mype,nprocs,mirror_blks,ix,iy,iz)



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

      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton
      use paramesh_interfaces, only     : morton_sort
      use paramesh_mpi_interfaces, only : morton_number

      implicit none
      include 'mpif.h'


      integer, intent(in)    ::  mype,nprocs,ix,iy,iz
      integer, intent(inout)   ::  mirror_blks(:,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables

      integer :: lb,mb,i,j,k,kk
      integer :: mort(6),istack


!??? Is 1000 the best value ?
      integer :: neigh_morts(6,2,npts_neigh)

      integer :: t1neigh_morts(6,npts_neigh)
      integer :: t2neigh_morts(6,npts_neigh)
      integer :: indx(npts_neigh),itemp
      integer :: inv_indx(npts_neigh)

      integer :: i_pe,j_pe,rem_block,rem_pe
      integer :: no_of_remote_neighs

      real    :: xmin,xmax,ymin,ymax,zmin,zmax
      real    :: coordm(3)
      logical :: lswap
      logical :: morton_less_than
      logical :: morton_greater_than
      logical :: morton_equal


      integer :: no_of_comms_to_send,kstack,max_no_to_be_received
      integer :: iprocs,isize,isrc,idest,itag,tag_offset,ll
      integer :: ierr,ierror,ierrorcode
      integer :: allocation_status
      integer :: max_no_of_blocks
      integer,dimension (:),  allocatable :: recvrequest
      integer,dimension (:,:),allocatable :: recvstatus


      common/track/ icomm_indicat
      integer :: icomm_indicat

      real    :: x0,y0,z0,bbsize(3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!
!
! This routine assumes that the grid blocks are ordered by morton
! number and that any blocks with different refinement levels but
! the same morton number are ordered from coarse to fine.


!--------------------------------------------------

      if(ix.le.0.and.iy.le.0.and.iz.le.0.and.mype.eq.0) then
        write(*,*)  & 
     &   'Paramesh error : arguments for amr_mirror_blks are wrong!', & 
     &   ' you must set one of ix,iy or iz to a value of 1 or 2.'
        call amr_abort()
      endif
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
      no_of_remote_neighs = 0

#ifdef DEBUG
      write(*,*) 'lrefine_max ',lrefine_max
      do lb=1,lnblocks
      if(nodetype(lb).eq.1) then
          x0 = coord(1,lb) -xmin
          y0 = coord(2,lb) -ymin
          z0 = coord(3,lb) -zmin
          bbsize(:) = bsize(:,lb)
          call morton_number( & 
     &          x0,y0,z0,bbsize,ndim,lrefine_max,lrefine(lb), & 
     &          mort)
          write(*,*) 'pe ',mype,' blk ',lb,' mortonnumber ', & 
     &          mort,' xyz0 ',x0,y0,z0

      endif
      enddo
#endif /* DEBUG */
!--------------------------------------------------
!
! Construct a list of morton numbers of any off-processor 
! mirror blocks.

      istack = 0

#ifdef DEBUGX
      write(*,*) 'xmin,ymin,zmin,xmax,ymax,zmax ', & 
     & xmin,ymin,zmin,xmax,ymax,zmax
#endif /* DEBUG */

!-------------
      do lb=1,lnblocks
      if(nodetype(lb).eq.1) then


! Mirror coordinate location offset from (xmin,ymin,zmin)

! equatorial symmetry
      if(ix.ne.0) then
        coordm(1) = grid_xmax - (coord(1,lb) - grid_xmin)
      else
        coordm(1) = coord(1,lb)
      endif
      if(iy.ne.0) then
        coordm(2) = grid_ymax - (coord(2,lb) - grid_ymin)
      else
        coordm(2) = coord(2,lb)
      endif
      if(iz.ne.0) then
        coordm(3) = grid_zmax - (coord(3,lb) - grid_zmin)
      else
        coordm(3) = coord(3,lb)
      endif


      coordm(1) = coordm(1)  - xmin
      if(ndim.ge.2) coordm(2) = coordm(2)  - ymin
      if(ndim.eq.3) coordm(3) = coordm(3)  - zmin

#ifdef DEBUG
      write(*,*) 'mirror_blks: pe ',mype,' lb ',lb, & 
     &    ' coord ',coord(:,lb),' mirror coordm ',coordm
#endif /* DEBUG */
      
! compute morton no for mirror block
      call morton_number(coordm(1),coordm(2),coordm(3), & 
     &                   bsize(1,lb),ndim, & 
     &                   lrefine_max,lrefine(lb), & 
     &                   mort)
      if( ( morton_less_than(mort(1:6), & 
     &      morton_limits(1:6,1,1,mype+1))  .or. & 
     &    ( morton_equal(mort(1:6), & 
     &      morton_limits(1:6,1,1,mype+1)) .and. & 
     &     (lrefine(lb).lt.morton_limits(6,2,1,mype+1)) ) ) & 
     &                    .or. & 
     &    ( morton_greater_than(mort(1:6), & 
     &                          morton_limits(1:6,1,2,mype+1))  .or. & 
     &    ( morton_equal(mort(1:6), & 
     &                    morton_limits(1:6,1,2,mype+1)) .and. & 
     &      (lrefine(lb).gt.morton_limits(6,2,2,mype+1)) ) ) & 
     &  ) then
        istack = istack + 1
        neigh_morts(:,1,istack) = mort(:)
        neigh_morts(:,2,istack) = lrefine(lb)
        mirror_blks(1,lb) = istack
        mirror_blks(2,lb) = -1000
#ifdef DEBUG
        write(*,*) 'pe ',mype,lb, & 
     &                  ' Mirror blk located off proc ', & 
     &     ' with mort no ',mort,' refine level ',lrefine(lb)
#endif /* DEBUG */
      else
        do mb = 1,lnblocks
          if( morton_equal(mort(1:6),mortonbnd(1:6,1,mb)).and. & 
     &        lrefine(lb).eq.mortonbnd(6,2,mb) ) & 
     &     then
            mirror_blks(1,lb) = mb
            mirror_blks(2,lb) = mype
#ifdef DEBUG
        write(*,*) 'pe ',mype,lb, & 
     &       'Mirror blk located on proc at ',mb
#endif /* DEBUG */
          endif
        enddo
      endif
      endif
      enddo

      no_of_remote_neighs = istack


!-------------
      if(no_of_remote_neighs.gt.0) then

!
! Sort this local list by morton no and refinement level


!
! Sort the list in order of increasing morton number, and then
! sort each sub-list with the same morton number in order of
! increasing refinement level.

! set indx to enable sort permutation to be stored in case we need
! it again.
       do i=1,istack
         indx(i) = i
       enddo

! sort the neigh_morts list according to their morton numbers
       do i = 1,6
          t1neigh_morts(i,1:istack) = neigh_morts(i,1,1:istack)
       end do
       call morton_sort(t1neigh_morts(:,1:istack),indx(1:istack),istack)
! now reorder the morton number part of neigh_morts
       do i = 1,6
          neigh_morts(i,1,1:istack) = t1neigh_morts(i,1:istack)
       end do
       t2neigh_morts(6,1:istack) = neigh_morts(6,2,1:istack)
       do i=1,istack
         neigh_morts(6,2,i) = t2neigh_morts(6,indx(i))
       enddo

! now update pointers to mirror block info after sorting.
!       do i=1,lnblocks
       do i=1,istack
         mirror_blks(1,indx(i)) = i
       enddo

#ifdef DEBUGX
       write(*,*) 'mirror pe ',mype,' after sorting neigh_morts ', & 
     &         neigh_morts(:,:,:istack)
#endif /* DEBUG  */

! now scan the list identifying segments with the same morton number
! and sort each segment based on refinement level
! order segments with same morton number in order of increasing
! refinement level
      lswap = .false.
      do while (lswap)
        lswap = .false.
        do i = 1,istack-1
          if(morton_equal(neigh_morts(1:6,1,i),neigh_morts(1:6,1,i+1)) & 
     &                        .and. & 
     &       neigh_morts(6,2,i).gt.neigh_morts(6,2,i+1) ) then
            lswap = .true.
            itemp = neigh_morts(6,2,i)
            neigh_morts(6,2,i) = neigh_morts(6,2,i+1)
            neigh_morts(6,2,i+1) = itemp
#ifdef DEBUGX
            write(*,*) 'pe ',mype,' swapping ',i,i+1
#endif /* DEBUG  */
          endif
        enddo
      enddo                           ! end do while


#ifdef DEBUGX
       write(*,*) 'mirror pe ',mype,' after sorting neigh_morts ', & 
     &         neigh_morts(:,:,:istack)
#endif /* DEBUG  */



!----------------------------------------------

!
! Initializations
      commatrix_send = 0
      commatrix_recv = 0
      pe_source = -1

!--------------------------------------------------
!
! Step 5.
! Construct a list of all processors from which the local processor should
! request morton number information.


! non-zero elements of COMMATRIX define which processor pairs need to 
! exchange morton number lists.

        commatrix_recv = 0
        do i = 1,no_of_remote_neighs
          i_pe = 1
          j_pe = -1
          do while(  & 
     &       (  morton_greater_than(neigh_morts(1:6,1,i), & 
     &                              morton_limits(1:6,1,2,i_pe)) & 
     &                               .or. & 
     &         (morton_equal(neigh_morts(1:6,1,i), & 
     &                       morton_limits(1:6,1,2,i_pe)).and. & 
     &          neigh_morts(6,2,i).gt.morton_limits(6,2,2,i_pe)  )  ) & 
     &          .and. (i_pe.le.nprocs) & 
     &            )
             i_pe = i_pe + 1
          enddo
          if(i_pe.le.nprocs) j_pe = i_pe
!
! If block has been located then update commatrix
          if(j_pe.ne.-1)  & 
     &      commatrix_recv(j_pe) =  commatrix_recv(j_pe) + 1

        enddo

#ifdef DEBUGY
        write(*,*) 'pe ',mype,' mirror commatrix bef gather ', & 
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
#ifdef DEBUGY
       write(*,*) 'pe ',mype,' no_of_comms_to_send ', & 
     &           no_of_comms_to_send
#endif /* DEBUG  */

!--------------------------------------------------
       endif                     ! end of no_of_remote_neighs if test
!--------------------------------------------------
       if(nprocs.gt.1) then
!--------------------------------------------------
!
! Step 6.
! provide the complete COMMATRIX to all processors

      call MPI_AlltoAll (commatrix_recv,       1,MPI_INTEGER, & 
     &                   commatrix_send,       1,MPI_INTEGER, & 
     &                   amr_mpi_meshComm,ierror)

#ifdef DEBUGX
        write(*,*) 'pe ',mype,' commatrix_send ', & 
     &             commatrix_send(1:nprocs)
 
      if(mype.eq.0) then
         write(*,'(" ")')
         write(*,'(" COMMUNICATION MATRIX: mirror")')
         write(*,'(" ")')
         write(*,'(" ",8i3)') (commatrix_send(i),i=1,nprocs)
         write(*,'(" ")')
      endif

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
       call MPI_ALLREDUCE(iprocs,  & 
     &                    max_no_to_be_received, & 
     &                    1, & 
     &                    MPI_INTEGER, & 
     &                    MPI_MAX, & 
     &                    amr_mpi_meshComm, & 
     &                    ierror)

#ifdef DEBUGX
       write(*,*) 'pe ',mype,' max_no_to_be_received ', & 
     &           max_no_to_be_received
#endif /* DEBUG  */


!--------------------------------------------------
!
! Step 8.
! Dynamically allocate memory to store the remote morton information.

       call MPI_ALLREDUCE(lnblocks,  & 
     &                    max_no_of_blocks, & 
     &                    1, & 
     &                    MPI_INTEGER, & 
     &                    MPI_MAX, & 
     &                    amr_mpi_meshComm, & 
     &                    ierror)

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

      tag_offset = 0

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
#ifdef DEBUGX
      write(*,*) 'pe ',mype,' mortonbnd receive of size ',isize, & 
     &          ' from ',isrc,' tag ',itag,' posted' & 
     &          ,' nprocs ',nprocs,' tag_offset ',tag_offset
#endif /* DEBUG  */
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
#ifdef DEBUGX
      write(*,*) 'pe ',mype,' mortonbnd of size ',isize, & 
     &          ' posting send to ',idest,' tag ',itag & 
     &          ,' nprocs ',nprocs,' tag_offset ',tag_offset
#endif /* DEBUG  */
             call MPI_Ssend(mortonbnd(1,1,1),isize,MPI_INTEGER, & 
     &            idest,itag,amr_mpi_meshComm,ierr)
#ifdef DEBUGX
      write(*,*) 'pe ',mype,' mortonbnd of size ',isize, & 
     &          ' sent to ',idest,' tag ',itag
#endif /* DEBUG  */
          endif
      enddo

      no_of_mortonbnds_received = k

#ifdef PM_UNIQUE_MPI_TAGS
      tag_offset = (nprocs-1)*nprocs + nprocs + tag_offset
#endif

      if(k.gt.0) & 
     &    call MPI_Waitall(k,recvrequest,recvstatus, & 
     &                     ierror)


#ifdef DEBUGX
      write(*,*) 'pe ',mype,' no_of_mortonbnds_received ', & 
     &          no_of_mortonbnds_received
#endif /* DEBUG  */
        
        
!--------------------------------------------------
!
! Loop over this processor^s list of required neighbor blocks,
! identifying their remote location from the morton information received

#ifdef DEBUGX 
          write(*,*) 'pe ',mype,' morton_limits ',morton_limits
#endif /* DEBUG  */
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
          enddo
          if(i_pe.le.nprocs) j_pe = i_pe

          rem_block = -1
          rem_pe = j_pe

          kk = -1
          do k=1,no_of_mortonbnds_received
            if(pe_source(k).eq.rem_pe) kk = k 
          enddo
#ifdef DEBUGX
          write(*,*) 'pe ',mype,' neigh i ',i,' i_pe ', & 
     &            i_pe,' neigh_morts ',neigh_morts(:,:,i), & 
     &             ' kk ',kk,' no_of_mortonbnds_received ', & 
     &              no_of_mortonbnds_received, & 
     &             ' rem_pe ',rem_pe, & 
     &             ' icomm_indicat ',icomm_indicat
#endif /* DEBUG  */

          if(kk.gt.0) then
!          do j=1,maxblocks
          do j=1,max_no_of_blocks
            if( morton_equal(r_mortonbnd(1:6,1,j,kk), & 
     &                       neigh_morts(1:6,1,i))  & 
     &                    .and. & 
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

          neigh_morts(6,1,i) = rem_block
          neigh_morts(6,2,i) = rem_pe-1

!#ifdef DEBUGX 
          write(*,*) 'pe ',mype,' neigh i ',i,' address ', & 
     &            neigh_morts(:,:,i)
!#endif /* DEBUG  */
        enddo

!--------------------------------------------------

      do i=1,istack
        inv_indx(indx(i)) = i
      enddo
      istack = 0

      do lb = 1,lnblocks
      if( (nodetype(lb).eq.1) .and. & 
     &    (mirror_blks(2,lb).eq.-1000) ) then

        istack = istack+1
        mirror_blks(1:2,lb) = neigh_morts(6,1:2,inv_indx(istack))
      else
        mirror_blks(1:2,lb) = -1
      endif
      enddo

!#ifdef DEBUGX
      do lb = 1,lnblocks
      if(nodetype(lb).eq.1) then
          write(*,*) 'pe ',mype,' lb ',lb,' mirror_blk ', & 
     &            mirror_blks(:,lb)
      endif
      enddo
      Call MPI_BARRIER(amr_mpi_meshComm, ierr)
!#endif /* DEBUG  */

!--------------------------------------------------
      endif                              ! end of nprocs > 1 if test
!--------------------------------------------------
!
! Error trapping.
      do lb = 1,lnblocks
      if(nodetype(lb).eq.1) then
      if(mirror_blks(1,lb).le.0.or.mirror_blks(2,lb).lt.0) then
          write(*,*) 'WARNING : Your grid is not symmetric. '
      endif
      endif
      enddo

      write(*,*) 'exiting mirror_blks'
      return
      end subroutine amr_mirror_blks
