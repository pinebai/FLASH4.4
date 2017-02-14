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
!#define DEBUGZ

      subroutine amr_refine_blocks (nprocs,mype)


!!****f* mpi_source/amr_refine_blocks
!! NAME
!!
!!   amr_refine_blocks
!! 
!! SYNOPSIS
!!
!!   call amr_refine_blocks (nprocs, mype)
!!
!!   call amr_refine_blocks (integer, integer)
!!
!! ARGUMENTS      
!!
!!   integer, intent(in) :: nprocs
!!     The number of processors being used.
!!
!!   integer, intent(in) :: mype     
!!     The calling processor number.
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
!!   paramesh_interfaces
!!   mpi_morton
!!   timings
!!   constants
!!   paramesh_mpi_interfaces
!!
!! CALLS
!!
!!   morton_neighbors
!!
!! RETURNS
!!
!!   This routine does not return anything.  Upon return new child blocks have been
!!   created and their locations in the tree data structure have been established.
!!
!! DESCRIPTION
!!
!!   The routine takes care of adding blocks to local list of blocks based 
!!   the setting of the flags set in the logical array 'refine(:)'.  If 
!!   refine(lb) is set to .true. for block lb, then the children of that block
!!   are created and added to the end of the list of on-processor blocks.
!!   The tree data is also updated by this routine to reflect these new blocks.
!!   The routine does NOT redistribute any data between processors.  That is handled
!!   by the routine amr_redist_blk.
!!
!!   Note that this routine is called by the routine 'amr_refine_derefine' and should
!!   not need to be called directly by a user's application.
!!
!! AUTHORS
!!
!!   Kevin Olson (original version 1996, MPI version 1997,1998)
!!   Modifications by Mike Zingale and Jonathan Dursi (U. of Chicago), 1999
!!
!!***

      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton
      use timings
      use constants
      use paramesh_mpi_interfaces, only :  morton_neighbors
      use paramesh_comm_data

      implicit none

      include 'mpif.h'

      integer, intent(in)    :: nprocs,mype

! local variables -------------------------------------------------------

      real :: hx,hy,hz
      real :: hh(mdim)      

      integer :: i,j,k,jj,ipar,ipar_proc,iproc
      integer :: ineigh,ineigh_proc
      integer :: lnblocks2,ichi
      integer :: ierr,jr,n3
      integer :: nsend,nrecv
      integer :: istart,iend,istart2,iend2,istep
      integer,allocatable :: nsend_to_proc(:), nrecv_pack(:)
      integer :: errorcode
      integer :: reqr(maxblocks_tr)
      integer :: kk(maxblocks_tr)
      integer :: nodetype_chi(nchild,maxblocks_tr)
      integer :: statr(MPI_STATUS_SIZE,maxblocks_tr)
      integer :: tneigh(2,mfaces,maxblocks_tr)
      integer :: idummy_array(1)
      integer :: iii,ii,jjj,kkk

      integer :: jr0,kr,isweep,nsweep
      logical :: lrecv
      logical :: lsend
      real    :: eps,accuracy

      logical :: refine_neigh(maxblocks_tr)
      double precision :: time1

! ------------------------------------------------------------------------
#ifdef DEBUG_FLOW_TRACE
       write(*,*) 'entered amr_refine_blocks'
#endif /* DEBUG_FLOW_TRACE */

      accuracy = 10./10.**precision(accuracy)
      eps = accuracy

#ifdef SAVE_MORTS
      real   :: xmin,ymin,zmin,xmax,ymax,zmax
      integer:: mort_neigh(6,3,3,3)

! Find the coordinate ranges (TO DO : This particular call need only be made at the start of the run,
! since it is always updated at the end of refine_derefine.)
      if (timing_mpi) then
         time1 = mpi_wtime()
      endif
      call mpi_amr_global_domain_limits
      if (timing_mpi) then
      timer_amr_global_domain_limits = & 
     &     timer_amr_global_domain_limits + mpi_wtime() - time1
      endif

! Compute xmin,ymin,zmin,xmax,ymax,zmax or get them from storage
      xmin = grid_xmin
      ymin = grid_ymin
      zmin = grid_zmin
      xmax = grid_xmax
      ymax = grid_ymax
      zmax = grid_zmax
#endif /* SAVE_MORTS */

#ifdef DEBUGX
      Call MPI_BARRIER(amr_mpi_meshComm, ierr)
      write(*,*) 'amr_refine_blocks-:-proc ',mype,' step 1'
      write(*,*) 'amr_refine_blocks : proc ',mype,' step 1'
      write(*,*) 'amr_refine_blocks : proc ',mype,' step 1'
      write(*,*) 'amr_refine_blocks : proc ',mype,' step 1'
      write(*,*) 'amr_refine_blocks : proc ',mype,' step 1'
      Call MPI_BARRIER(amr_mpi_meshComm, ierr)
      call amr_test_neigh_values()
#endif /* DEBUG */

! 1) refine blocks, create their children, the children's parents,
!    turn the children on, and turn the parents off

      do i = 1,maxblocks_tr
         newchild(i) = .FALSE.
      end do

!!!!!!!!!!!!!!

      lnblocks2 = lnblocks

      do i = 1,lnblocks

         if (refine(i)) then ! refine block 'i' on this processor

            do j = 1,nchild ! create 'nchild' new child blocks for block 'i'

               lnblocks2 = lnblocks2 + 1

      if(lnblocks2.gt.maxblocks_tr) then
       open(unit=30,file='amr_log',status='unknown', & 
     &        position='append')
       write(30,*) 'PARAMESH ERROR !'
       write(30,*) 'Too many blocks created! '
       write(30,*) 'Increase MAXBLOCKS_TR and rerun! '
       close(30)
       write(*,*) 'PARAMESH ERROR !'
       write(*,*) 'Too many blocks created! '
       write(*,*) 'Increase MAXBLOCKS_TR and rerun! '
       call MPI_ABORT(amr_mpi_meshComm,errorcode,ierr)
      endif

               child(1,j,i) = lnblocks2 ! child j's on-processor id
               child(2,j,i) = mype      ! child j's processor no.

               lrefine(lnblocks2) = lrefine(i) + 1 ! refinement level of child
               newchild(lnblocks2) = .TRUE. ! this is a new child

               parent(1,lnblocks2) = i ! child j's parent
               parent(2,lnblocks2) = mype ! child j's parent's processor

               bnd_box(:,:,lnblocks2) = bnd_box(:,:,i)

               ! ordering of children in space is morton ordered !

               hh(:) = bsize(:,i)/4.
               if (j.eq.1) then
                  hx = -1.
                  hy = -1.
                  hz = -1.
                  bnd_box(2,1,lnblocks2) = coord(1,i)
                  bnd_box(2,2,lnblocks2) = coord(2,i)
                  bnd_box(2,3,lnblocks2) = coord(3,i)
               else if (j.eq.2) then
                  hx = 1.
                  hy = -1.
                  hz = -1.
                  bnd_box(1,1,lnblocks2) = coord(1,i)
                  bnd_box(2,2,lnblocks2) = coord(2,i)
                  bnd_box(2,3,lnblocks2) = coord(3,i)
               else if (j.eq.3) then
                  hx = -1.
                  hy = 1.
                  hz = -1.
                  bnd_box(2,1,lnblocks2) = coord(1,i)
                  bnd_box(1,2,lnblocks2) = coord(2,i)
                  bnd_box(2,3,lnblocks2) = coord(3,i)
               else if (j.eq.4) then
                  hx = 1.
                  hy = 1.
                  hz = -1.
                  bnd_box(1,1,lnblocks2) = coord(1,i)
                  bnd_box(1,2,lnblocks2) = coord(2,i)
                  bnd_box(2,3,lnblocks2) = coord(3,i)
               else if (j.eq.5) then
                  hx = -1.
                  hy = -1.
                  hz = 1.
                  bnd_box(2,1,lnblocks2) = coord(1,i)
                  bnd_box(2,2,lnblocks2) = coord(2,i)
                  bnd_box(1,3,lnblocks2) = coord(3,i)
               else if (j.eq.6) then
                  hx = 1.
                  hy = -1.
                  hz = 1.
                  bnd_box(1,1,lnblocks2) = coord(1,i)
                  bnd_box(2,2,lnblocks2) = coord(2,i)
                  bnd_box(1,3,lnblocks2) = coord(3,i)
               else if (j.eq.7) then
                  hx = -1.
                  hy = 1.
                  hz = 1.
                  bnd_box(2,1,lnblocks2) = coord(1,i)
                  bnd_box(1,2,lnblocks2) = coord(2,i)
                  bnd_box(1,3,lnblocks2) = coord(3,i)
               else if (j.eq.8) then
                  hx = 1.
                  hy = 1.
                  hz = 1.
                  bnd_box(1,1,lnblocks2) = coord(1,i)
                  bnd_box(1,2,lnblocks2) = coord(2,i)
                  bnd_box(1,3,lnblocks2) = coord(3,i)
               end if

               if(ndim.lt.3) bnd_box(:,3,lnblocks2) = bnd_box(:,3,i)
               if(ndim.lt.2) bnd_box(:,2,lnblocks2) = bnd_box(:,2,i)

               which_child(lnblocks2) = j
                     
               coord(1,lnblocks2) = coord(1,i) + hx*hh(1)
               if (ndim.ge.2) then
                  coord(2,lnblocks2) = coord(2,i) + hy*hh(2)
               end if
               if (ndim.eq.3) then
                  coord(3,lnblocks2) = coord(3,i) + hz*hh(3)
               end if

               bsize(:,lnblocks2) = bsize(:,i)/2.

               bflags(:,lnblocks2) = bflags(:,i)

       if (empty_cells) then
          if(empty(i).eq.1) empty(lnblocks2)=1
       end if
#ifdef SAVE_MORTS
      call morton_neighbors(xmin,ymin,zmin,xmax,ymax,zmax, & 
     &                      lperiodicx,lperiodicy,lperiodicz, & 
     &                      coord(:,lnblocks2),bsize(:,lnblocks2),ndim, & 
     &                      lrefine(lnblocks2),lrefine_max, & 
     &                      mort_neigh, & 
     &                      bnd_box(1,1,lnblocks2))
      surr_morts(:,:,1:1+2*k2d,1:1+2*k3d,lnblocks2) =  & 
     &                           mort_neigh(:,:,2-k2d:2+k2d,2-k3d:2+k3d)
      write(*,*) 'pe ',mype,' added surr_morts for new block ', & 
     &  lnblocks2, & 
     & ' mort_neigh(6,:,2-k2d:2+k2d,2-k3d:2+k3d) ', & 
     &   mort_neigh(6,:,2-k2d:2+k2d,2-k3d:2+k3d)
#endif /* SAVE_MORTS */

            end do

       if (empty_cells) then
          if(empty(i).eq.1) empty(i)=0
       end if

         end if

       end do

! Connect neighbors of newly created sub-blocks
!
!                    4              6
!                    |              |
!  in x-y plane  1 - i - 2 ; in z   i
!                    |              |
!                    3              5
!

! In some cases, old data in the parts of the neigh array
! about to be used by the new children can confuse the algorithm
! which searches for off processor neighbors. In particular, if
! an external boundary is not properly detected, then a spurious
! neighbor link can be set up in section XXX below.
! Reinitializing these parts of neigh avoids the confusion.
      neigh(:,:,lnblocks+1:lnblocks2) = -1
      child(:,:,lnblocks+1:lnblocks2) = -1


! 1) connect with siblings (which at this point are all on processor)
      
      do i = 1,lnblocks
         if (refine(i)) then
            do j = 1,nchild ! cycle through children
               ichi = child(1,j,i)
               do k = 1,nfaces
                  neigh(2,k,ichi) = mype
               end do
            end do

            ! connect in x direction
               
            do j = 2,nchild,2
               ichi = child(1,j,i)
               k = j - 1
               ineigh = child(1,k,i)
               neigh(1,1,ichi) = ineigh
            end do
            do j = 1,nchild-1,2
               ichi = child(1,j,i)
               k = j + 1
               ineigh = child(1,k,i)
               neigh(1,2,ichi) = ineigh
            end do

            ! connect in y direction

            if(ndim.ge.2) then
            do j = 3,4
               ichi = child(1,j,i)
               k = j - 2
               ineigh = child(1,k,i)
               neigh(1,3,ichi) = ineigh
            end do
            if (ndim.eq.3) then
            do j = 7,8
               ichi = child(1,j,i)
               k = j - 2
               ineigh = child(1,k,i)
               neigh(1,3,ichi) = ineigh
            end do
            end if
            do j = 1,2
               ichi = child(1,j,i)
               k = j + 2
               ineigh = child(1,k,i)
               neigh(1,4,ichi) = ineigh
            end do
            if (ndim.eq.3) then
            do j = 5,6
               ichi = child(1,j,i)
               k = j + 2
               ineigh = child(1,k,i)
               neigh(1,4,ichi) = ineigh
            end do
            end if
            end if

            if (ndim == 3) then

            ! connect in z direction

            if (ndim.eq.3) then
            do j = 5,8
               ichi = child(1,j,i)
               k = j - 4
               ineigh = child(1,k,i)
               neigh(1,5,ichi) = ineigh
            end do
            do j = 1,4
               ichi = child(1,j,i)
               k = j + 4
               ineigh = child(1,k,i)
               neigh(1,6,ichi) = ineigh
            end do
            end if

            end if ! end if (ndim == 3)
         end if
      end do


! 2) connect with off-processor neighbors by looking at neighbors of parent

! In order to establish pointers between children of newly refined blocks,
! each parent block must know if its neighbor is refined. In the next section
! the refine flags are passed between neighbors. In a regular cartesian grid
! this is a simple send-recv operation across every face. For face 3 say
! (ie j=3) each block sends its refine value to its face 3 neighbor. Each
! block also issues a receive from its face 4 neighbor. In this way sends and receives
! are automatically matched.
! In non-cartesian geometry with singular lines, eg a pole in spherical coordinates,
! block faces touching the pole, touch their neighbor with the same face number,
! eg face 3 touches face 3. Therefore in the send/recv example above, polar blocks
! must issue a send but also a receive for their face 3 neighbor, not their face 4
! neighbor.
! To make this work, we have written the following section in two sweeps. The
! first sweep operates as in the cartesian example, except blocks with polar faces
! do not issue any sends or receives for their polar faces. On the second sweep,
! only the polar faces issue sends and receives.

      tneigh(:,:,:) = 0

! Send refine flags to neighbors

      do j = 1,nfaces

      nsweep = 1

      if(spherical_pm) then
         if(j.eq.3.or.j.eq.4) nsweep = 2
      endif

      do isweep = 1,nsweep

         if (j.eq.1) then
            jr = 2
            istart = 1
            iend = nchild-1
            istart2 = 2
            iend2 = nchild
            istep = 2
         elseif (j.eq.2) then
            jr = 1
            istart = 2
            iend = nchild
            istart2 = 1
            iend2 = nchild-1
            istep = 2
         elseif (j.eq.3) then
            jr = 4               
            istart = 1
            iend = 2
            istart2 = 3
            iend2 = 4
            istep = 1
         elseif (j.eq.4) then
            jr = 3
            istart = 3
            iend = 4
            istart2 = 1
            iend2 = 2
            istep = 1
         elseif (j.eq.5) then
            jr = 6
            istart = 1
            iend = 4
            istart2 = 5
            iend2 = 8
            istep = 1
         elseif (j.eq.6) then
            jr = 5
            istart = 5
            iend = 8
            istart2 = 1
            iend2 = 4
            istep = 1
         end if

         refine_neigh(:) = .FALSE.
         refine(lnblocks+1:maxblocks_tr) = .FALSE.

#ifdef DEBUG
        Call MPI_BARRIER(amr_mpi_meshComm, ierr)
        write(*,*) 'pe ',mype,' lnblocks lnblocks2 ',lnblocks,lnblocks2
        Call MPI_BARRIER(amr_mpi_meshComm, ierr)
#endif /* DEBUG */

         nrecv = 0
         do i = 1,lnblocks

! do receive from opposite index face, unless opposite face is
! a theta face in spherical coords
         lrecv=.true.

         if(isweep.eq.1) then

         if(spherical_pm) then
! In first sweep switch off any receives which might be generated at polar faces
          if(j.eq.3.and.abs(bnd_box(2,2,i)-pi).lt.eps) lrecv=.false.
          if(j.eq.4.and.abs(bnd_box(1,2,i)).lt.eps) lrecv=.false.
         endif

            if (lrecv) then

            if (neigh(1,jr,i).gt.0) then
               if (neigh(2,jr,i).ne.mype) then
                  nrecv = nrecv + 1
#ifdef DEBUG
      write(*,*) 'amr_refine_block : proc ',mype,' blk ',i, & 
     &           ' nrecv ',nrecv, & 
     &           ' source pe ',neigh(2,jr,i), & 
     &           ' posting ref normal recv tag ',neigh(1,jr,i),jr
#endif /* DEBUG */
                  call MPI_logical_IRECV(refine_neigh(i), & 
     &                           1, & 
     &                           MPI_LOGICAL, & 
     &                           neigh(2,jr,i), & 
     &                           neigh(1,jr,i),         & ! message tag ?
     &                           amr_mpi_meshComm, & 
     &                           reqr(nrecv), & 
     &                           ierr)
               end if
            end if

            end if

           endif               ! isweep = 1

         if(spherical_pm) then
         if(isweep.eq.2) then
! This section executes receives for polar faces only.
         lrecv = .false.
         if(j.eq.4.and.abs(bnd_box(2,2,i)-pi).lt.eps) lrecv=.true.
         if(j.eq.3.and.abs(bnd_box(1,2,i)).lt.eps) lrecv=.true.
         if(abs(bnd_box(1,2,i)).lt.eps.and. & 
     &        abs(bnd_box(2,2,i)-pi).lt.eps) then
             lrecv=.true.
          endif
            if (lrecv) then
            if (neigh(1,j,i).gt.0) then
               if (neigh(2,j,i).ne.mype) then
                  nrecv = nrecv + 1
                  call MPI_logical_IRECV(refine_neigh(i), & 
     &                           1, & 
     &                           MPI_LOGICAL, & 
     &                           neigh(2,j,i), & 
     &                           neigh(1,j,i),         & ! message tag ?
     &                           amr_mpi_meshComm, & 
     &                           reqr(nrecv), & 
     &                           ierr)
#ifdef DEBUG
      write(*,*) 'amr_refine_block : proc ',mype,' blk ',i, & 
     &           ' nrecv ',nrecv, & 
     &           ' source pe ',neigh(2,j,i), & 
     &           ' posting ref recv tag ',neigh(1,j,i),j
#endif /* DEBUG */
               end if
            end if

            end if
           endif               ! isweep = 2
         endif  ! end if (spherical_pm
         end do  
#ifdef DEBUG
      write(*,*) 'amr_refine_block : proc ',mype,' final nrecv ',nrecv
      call amr_flush(6)
      Call MPI_BARRIER(amr_mpi_meshComm, ierr)
#endif /* DEBUG */

         nsend = 0
         do i = 1,lnblocks

         lsend = .true.

         if(isweep.eq.1) then
         if(spherical_pm) then
! In first sweep switch off any sends which might be generated at polar faces
         if(j.eq.4.and.abs(bnd_box(2,2,i)-pi).lt.eps) lsend=.false.
         if(j.eq.3.and.abs(bnd_box(1,2,i)).lt.eps) lsend=.false.
         endif

            if (lsend) then
            ineigh = neigh(1,j,i)
            ineigh_proc = neigh(2,j,i)
            if (ineigh.ge.1) then
               if (ineigh_proc.ne.mype) then
                  nsend = nsend + 1
#ifdef DEBUG
      write(*,*) 'amr_refine_block : proc ',mype,' nsend ',nsend, & 
     &           ' dest pe ',ineigh_proc, & 
     &           ' posting ref send tag ',i,' j ',j,' ineigh ',ineigh
#endif /* DEBUG */
                  call MPI_logical_SSEND(refine(i), & 
     &                           1, & 
     &                           MPI_LOGICAL, & 
     &                           ineigh_proc, & 
     &                           i,                      & ! message tag ?
     &                           amr_mpi_meshComm, & 
     &                           ierr)
               else
                  if(lsend) refine_neigh(ineigh) = refine(i)
               end if
            end if
            end if
         endif               ! isweep = 1

         if(isweep.eq.2) then
         if(spherical_pm) then
         lsend = .false.
! In second sweep switch on sends associated with polar faces only
         if(j.eq.4.and.abs(bnd_box(2,2,i)-pi).lt.eps) lsend=.true.
         if(j.eq.3.and.abs(bnd_box(1,2,i)).lt.eps) lsend=.true.
         endif

            if (lsend) then
            ineigh = neigh(1,j,i)
            ineigh_proc = neigh(2,j,i)
            if (ineigh.ge.1) then
               if (ineigh_proc.ne.mype) then
                  nsend = nsend + 1
#ifdef DEBUG
      write(*,*) 'amr_refine_block : proc ',mype,' nsend ',nsend, & 
     &           ' dest pe ',ineigh_proc, & 
     &           ' posting ref send tag ',i,' j ',j,' ineigh ',ineigh
#endif /* DEBUG */
                  call MPI_logical_SSEND(refine(i), & 
     &                           1, & 
     &                           MPI_LOGICAL, & 
     &                           ineigh_proc, & 
     &                           i,                      & ! message tag ?
     &                           amr_mpi_meshComm, & 
     &                           ierr)
               else
                  if(lsend) refine_neigh(ineigh) = refine(i)
               end if
            end if
            end if
         endif               ! isweep = 2


         end do
#ifdef DEBUG
      write(*,*) 'amr_refine_block : proc ',mype,' final nsend ',nsend
      call amr_flush(6)
      Call MPI_BARRIER(amr_mpi_meshComm, ierr)
#endif /* DEBUG */

         if (nrecv.gt.0) then
            call MPI_WAITALL(nrecv,reqr,statr,ierr)
         end if
 100     continue

#ifdef DEBUG
      Call MPI_BARRIER(amr_mpi_meshComm, ierr)
      write(*,*) 'pe ',mype,' face ',j,' refine_neigh ', & 
     &                refine_neigh(1:lnblocks)
      write(*,*) 'pe-',mype,' at label 100 again'
      write(*,*) 'pe ',mype,' at label 100 again'
      write(*,*) 'pe ',mype,' at label 100 again'
      call amr_flush(6)
      Call MPI_BARRIER(amr_mpi_meshComm, ierr)
#endif /* DEBUG */

         do jj = istart,iend,istep
            
            if (j.eq.1) then
               k = jj+1         ! child no. of neighbor which lies on border
            elseif (j.eq.2) then
               k = jj-1
            elseif (j.eq.3) then
               k = jj+2
            elseif (j.eq.4) then
               k = jj-2
            elseif (j.eq.5) then
               k = jj+4
            elseif (j.eq.6) then
               k = jj-4
            end if

         nrecv = 0
         do i = 1,lnblocks
         kr = k
            if (refine(i)) then

              if(isweep.eq.1) then
              lrecv = .true.

      if(spherical_pm) then
! In first sweep switch off any receives which might be generated at polar faces
         if(j.eq.3.and.abs(bnd_box(1,2,i)).lt.eps) lrecv=.false.
         if(j.eq.4.and.abs(bnd_box(2,2,i)-pi).lt.eps) lrecv=.false.
      endif

              elseif(isweep.eq.2) then
              lrecv = .false.

      if(spherical_pm) then
! In first sweep switch off any receives which might be generated at polar faces
         if(j.eq.3.and.abs(bnd_box(1,2,i)).lt.eps) lrecv=.true.
         if(j.eq.4.and.abs(bnd_box(2,2,i)-pi).lt.eps) lrecv=.true.
      endif

              endif

              if (lrecv) then

               ineigh = neigh(1,j,i) ! neighbor of parent
               ineigh_proc = neigh(2,j,i) ! neighbor of parent
               if (ineigh.gt.0) then
                                            ! 'i' on left side
                  ichi = child(1,jj,i) ! child 'jj' of parent 'i'

                     ! fetch child no. k of ineigh from ineigh_proc
                     ! and store in neighbor 1 of child jj (ichi)

                  if (ineigh_proc.ne.mype) then
                     nrecv = nrecv + 1
                     call MPI_int_IRECV(neigh(1,j,ichi), & 
     &                              2, & 
     &                              MPI_INTEGER, & 
     &                              neigh(2,j,i), & 
     &                              neigh(1,j,i), & 
     &                              amr_mpi_meshComm, & 
     &                              reqr(nrecv), & 
     &                              ierr)
#ifdef DEBUG
      write(*,*) 'amr_refine_block : sweep ',isweep, & 
     &           ' proc ',mype,' blk ',i, & 
     &           ' jface ',j, & 
     &           ' nrecv ',nrecv, & 
     &           ' source pe ',neigh(2,j,i), & 
     &           ' posting child recv sourceblk/tag ',neigh(1,j,i)
#endif /* DEBUG */
                  else

      if(spherical_pm) then
                     if(isweep.eq.2) then
! if this child is at the pole then modify the index of the child
! which should send back.
! (ie if spherical pole block then k=jj in next line)
                     if( j.eq.3 & 
     &                   .and.abs(bnd_box(1,2,i)).lt.eps) then
                       kr = jj
                     endif
                     if( j.eq.4 & 
     &                   .and.abs(bnd_box(2,2,i)-pi).lt.eps) then
                       kr = jj
                     endif
                     endif                 ! isweep.eq.2
      endif

                     neigh(:,j,ichi) = child(:,kr,ineigh)
                  end if
               end if
              end if
            end if
         end do

#ifdef DEBUG
      write(*,*) 'amr_refine_block : proc ',mype,' face ',j, & 
     &           ' child final nrecv ',nrecv
      Call MPI_BARRIER(amr_mpi_meshComm, ierr)
#endif /* DEBUG */

         nsend = 0
         do i = 1,lnblocks
         kr = k

            if (refine_neigh(i)) then

         if(isweep.eq.1) then

      if(spherical_pm) then
         lsend = .true.
! In first sweep switch off sends associated with polar faces
         if(j.eq.4.and.abs(bnd_box(1,2,i)).lt.eps) lsend=.false.
         if(j.eq.3.and.abs(bnd_box(2,2,i)-pi).lt.eps) lsend=.false.
      endif

            if (lsend) then

               ineigh = neigh(1,jr,i) ! right neighbor of parent
               ineigh_proc = neigh(2,jr,i)
               if (ineigh.gt.0) then
                  ! send child no. jj of block i to neigh on right
                  if (ineigh_proc.ne.mype) then
                     nsend = nsend + 1
#ifdef DEBUG
      write(*,*) 'amr_refine_block : sweep ',isweep, & 
     &           ' proc ',mype,' nsend ',nsend, & 
     &           ' dest pe ',neigh(2,jr,i), & 
     &           ' posting basic child send tag ',i
#endif /* DEBUG */
                     call MPI_int_SSEND(child(1,kr,i), & 
     &                    2, & 
     &                    MPI_INTEGER, & 
     &                    neigh(2,jr,i), & 
     &                    i, & 
     &                    amr_mpi_meshComm, & 
     &                    ierr)
                  end if
               end if

             end if
         endif


         if(isweep.eq.2) then

      if(spherical_pm) then
         lsend = .false.
! In second sweep switch on sends associated with polar faces only
         if(j.eq.3.and.abs(bnd_box(1,2,i)).lt.eps) lsend=.true.
         if(j.eq.4.and.abs(bnd_box(2,2,i)-pi).lt.eps) lsend=.true.
            if (lsend) then

               jr0 = j
               ineigh = neigh(1,jr0,i) ! right neighbor of parent
               ineigh_proc = neigh(2,jr0,i)
               if (ineigh.gt.0 .and. j.eq.jr0) then
                  ! send child no. jj of block i to neigh on right
                  if (ineigh_proc.ne.mype) then
                     nsend = nsend + 1

! if spherical pole block then kr=jj in next line
                     if( j.eq.3 & 
     &                   .and.abs(bnd_box(1,2,i)).lt.eps) then
                       kr = jj
                     endif
                     if( j.eq.4 & 
     &                   .and.abs(bnd_box(2,2,i)-pi).lt.eps) then
                       kr = jj
                     endif
#ifdef DEBUG
      write(*,*) 'amr_refine_block : sweep ',isweep, & 
     &           ' proc ',mype,' nsend ',nsend, & 
     &           ' dest pe ',neigh(2,jr0,i), & 
     &           ' dest blk ',neigh(1,jr0,i), & 
     &           ' posting polar child send tag ',i
#endif /* DEBUG */
                     call MPI_int_SSEND(child(1,kr,i), & 
     &                    2, & 
     &                    MPI_INTEGER, & 
     &                    neigh(2,jr0,i), & 
     &                    i, & 
     &                    amr_mpi_meshComm, & 
     &                    ierr)
                  end if
               end if
            end if
      endif

         endif


            end if
         end do
#ifdef DEBUG
      write(*,*) 'amr_refine_block : proc ',mype,' face ',j, & 
     &           ' final child nsend ',nsend
#endif /* DEBUG */
         
         if (nrecv.gt.0) then
            call MPI_WAITALL(nrecv,reqr,statr,ierr)
         end if


         end do ! end loop over jj

         if (j.eq.3.and.istart.eq.1.and.ndim.eq.3) then
            istart = 5
            iend = 6
            istart2 = 7
            iend2 = 8
            go to 100
         elseif (j.eq.4.and.istart.eq.3.and.ndim.eq.3) then
            istart = 7
            iend = 8
            istart2 = 5
            iend2 = 6
            go to 100
         end if

! 1) count no. of sends on each proc. to all other procs
         allocate(nsend_to_proc(0:nprocs))

         nsend_to_proc(:) = 0
         do i = lnblocks+1,lnblocks2
            if (newchild(i)) then
               ineigh = neigh(1,j,i)
               ineigh_proc = neigh(2,j,i)
               if (ineigh.ge.1) then
                  if (ineigh_proc.ne.mype) then
                     nsend_to_proc(ineigh_proc) =  & 
     &                    nsend_to_proc(ineigh_proc) + 1
                  end if
               end if
            end if
         end do

! 2) collect data for `this' proc from other procs so that
!    the total no. of receives to post can be computed

         allocate(nrecv_pack(0:nprocs-1))
! Change made by M. Zingale and J. Dursi -------------------------------------
         call MPI_int_AllReduce(nsend_to_proc, nrecv_pack, nprocs, & 
     &        MPI_INTEGER, MPI_SUM, amr_mpi_meshComm, ierr)
         nrecv = nrecv_pack(MyPE)
         deallocate(nrecv_pack)
         deallocate(nsend_to_proc)

!-----------------------------------------------------------------------------

         do i = 1,nrecv
#ifdef DEBUGZ
           write(*,*) 'DZ: pe ',mype,' posting recv'
#endif /* DEBUGZ */
            call MPI_int_IRECV(kk(i), & 
     &                     1, & 
     &                     MPI_INTEGER, & 
     &                     MPI_ANY_SOURCE, & 
     &                     MPI_ANY_TAG, & 
     &                     amr_mpi_meshComm, & 
     &                     reqr(i), & 
     &                     ierr)
         end do  

! section XXX
! NOW new children send to their new neighbors so that they can be set
                     
         nsend = 0
         do i = lnblocks+1,lnblocks2
            if (newchild(i)) then
               ineigh = neigh(1,j,i)
               ineigh_proc = neigh(2,j,i)
               if (ineigh.ge.1) then
                  if (ineigh_proc.ne.mype) then
                     nsend = nsend + 1
                     idummy_array(1) = ineigh
#ifdef DEBUGZ
           write(*,*) 'DZ: pe ',mype,' posting send to ', & 
     &                  ineigh_proc,' tag ',i
#endif /* DEBUGZ */
                     call MPI_int_SSEND(idummy_array(1), & 
     &                              1, & 
     &                              MPI_INTEGER, & 
     &                              ineigh_proc, & 
     &                              i, & 
     &                              amr_mpi_meshComm, & 
     &                              ierr)
                  else
                     tneigh(1,jr,ineigh) = i
                     tneigh(2,jr,ineigh) = mype
                  end if
               end if
            end if
         end do

! it seems to work only if tneigh is used ????

         if (nrecv.gt.0) then
            call MPI_WAITALL(nrecv,reqr,statr,ierr)
            do i = 1,nrecv
               tneigh(1,jr,kk(i)) = statr(MPI_TAG,i)
               tneigh(2,jr,kk(i)) = statr(MPI_SOURCE,i)
            end do
         end if
         call MPI_BARRIER(amr_mpi_meshComm,ierr) ! NEEDED ????

      enddo                  ! end of isweep loop

      end do ! end loop over faces

                     
      do i = 1,lnblocks2
         do j = 1,nfaces
            if (tneigh(1,j,i).gt.0.and..not.newchild(i).and. & 
     &           tneigh(1,j,i).ne.neigh(1,j,i)) then
               neigh(1,j,i) = tneigh(1,j,i)
               neigh(2,j,i) = tneigh(2,j,i)
          
            end if
         end do
      end do
      

      lnblocks = lnblocks2

! reset node types
      
      do i = 1,lnblocks
         if (newchild(i)) nodetype(i) = 1
         if (refine(i)) then
            nodetype(i) = 2
         end if
      end do

      nrecv = 0
      nodetype_chi(:,1:lnblocks) = 0
      do i = 1,lnblocks
         do j = 1,nchild
            if (child(1,j,i).gt.0) then
            if (child(2,j,i).ne.mype) then
               nrecv = nrecv + 1
               call MPI_int_IRECV(nodetype_chi(j,i), & 
     &                        1, & 
     &                        MPI_INTEGER, & 
     &                        child(2,j,i), & 
     &                        child(1,j,i), & 
     &                        amr_mpi_meshComm, & 
     &                        reqr(nrecv), & 
     &                        ierr)
            else
               nodetype_chi(j,i) = nodetype(child(1,j,i))
            end if
            end if
         end do
      end do  

      nsend = 0
      do i = 1,lnblocks
         if (parent(1,i).gt.0) then
         if (parent(2,i).ne.mype) then
            nsend = nsend + 1
            call MPI_int_SSEND(nodetype(i), & 
     &                     1, & 
     &                     MPI_INTEGER, & 
     &                     parent(2,i), & 
     &                     i, & 
     &                     amr_mpi_meshComm, & 
     &                     ierr)
         end if
         end if
      end do

      if (nrecv.gt.0) then
         call MPI_WAITALL(nrecv,reqr,statr,ierr)
      end if
      do i = 1,lnblocks
         n3 = 0
         do j = 1,nchild
            if (nodetype_chi(j,i).ne.1) then
               n3 = n3 + 1
            end if
            if (n3.eq.nchild) nodetype(i) = 3
         end do
      end do 

! Now set neighbor pointers of new children if they lie on a boundary

      do i = 1,lnblocks

         if (newchild(i)) then

            ! fetch i's parent

            ipar = parent(1,i)
            ipar_proc = parent(2,i)

            do j = 1,nfaces

               if (neigh(1,j,i).le.-1) then ! this neighbor may be on a border

                  ! fetch i's parent's neighbor j

! Here we know that all new children on the same processor as
! the parent so we don't need any communications !!!

                  ineigh = neigh(1,j,ipar)
               
                  ! if parent's neighbor is lt -1 then i's neighbor is
                  ! also on the domain border and is set to the parent's
                  ! value

                  if (ineigh.le.-20) neigh(1,j,i) = ineigh

               end if

            end do

         end if

      end do

! reset refine flags

      do i = 1,maxblocks_tr
         refine(i) = .FALSE.
      end do

#ifdef DEBUG
      Call MPI_BARRIER(amr_mpi_meshComm, ierr)
      write(*,*) 'pe-',mype,'-near end of amr_refine_blocks'
      write(*,*) 'pe-',mype,'-calling amr_test_neigh_values'
!      call amr_test_neigh_values()
      do iproc=0,nprocs-1
      Call MPI_BARRIER(amr_mpi_meshComm, ierr)
      if(iproc.eq.mype) then
      do i=1,lnblocks
      write(*,*) 'proc ',mype,' neigh ',i,neigh(:,1:2*ndim,i)
      enddo
      endif
      enddo
      write(*,*) 'pe-',mype,'-exited amr_test_neigh_values'
      write(*,*) 'pe ',mype,' exiting amr_refine_blocks'
      Call MPI_BARRIER(amr_mpi_meshComm, ierr)
#endif /* DEBUG */
#ifdef DEBUG_FLOW_TRACE
       write(*,*) 'exiting amr_refine_blocks'
#endif /* DEBUG_FLOW_TRACE */

      return
      end subroutine amr_refine_blocks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine amr_check_refine(nprocs,mype,icontinue)

!!****f* mpi_source/amr_check_refine
!! NAME
!!
!!   amr_check_refine
!! 
!! SYNOPSIS
!!
!!   call amr_refine_blocks (nprocs, mype, icontinue)
!!
!!   call amr_refine_blocks (integer, integer, integer)
!!
!! ARGUMENTS      
!!
!!   integer, intent(in) :: nprocs
!!     The number of processors being used.
!!
!!   integer, intent(in) :: mype     
!!     The calling processor number.
!!
!!   integer, intent(in) :: icontinue     
!!     A variable controling algorithm behavior.
!!
!! INCLUDES
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   constants
!!
!! CALLS
!!
!!   morton_neighbors
!!
!! RETURNS
!!
!!   This routine does not return anything.  Upon return additional refine flags
!!   may be set to .true. to ensure a valid mesh.
!!
!! DESCRIPTION
!!
!!   This routine performs checks on the list of refine flags to make sure that
!!   a valid mesh results when the new, child blocks are added.  It makes sure that
!!   a jump in refinement of no more than one level is ensured.
!!  
!!   The routine is called only my amr_refine_derefine and a user's application should
!!   not need to call this routine directly.
!!
!! AUTHORS
!!
!!   Kevin Olson (1999)
!!
!!***

      use paramesh_dimensions
      use physicaldata
      use tree
      use constants
      use paramesh_comm_data

      implicit none

      integer, intent(in)    :: nprocs,mype
      integer, intent(in)    :: icontinue

      include 'mpif.h'

! local variables ----------------------------------------------------------

      integer :: i,j,ipar,ipar_proc
      integer :: ineigh,ineigh_proc
      integer :: ix, iy, iz, ix_n, iy_n, iz_n
      integer :: ierr,jr
      integer :: nsend,nrecv
      integer :: reqr(maxblocks_tr)
      integer :: statr(MPI_STATUS_SIZE,maxblocks_tr)

      logical :: repeat,repeat_t
      logical :: ref_testt(2,1+k2d,1+k3d,maxblocks_tr)
      logical :: ref_testts(2,1+k2d,1+k3d,maxblocks_tr)
      logical :: ref_test(2,1+k2d,1+k3d,maxblocks_tr)
      logical :: lreduce_datain(1), lreduce_dataout(1)

      real, parameter :: eps = 1.e-10
      integer :: jr0
      logical :: lrecv,lsend

! --------------------------------------------------------------------------
#ifdef DEBUG_FLOW_TRACE
         write(*,*) 'pe ',mype,' entered amr_check_refine'
#endif /* DEBUG_FLOW_TRACE */


! Turn off refine flags of non-leaf blocks
      do i = 1,lnblocks
         if (nodetype(i).gt.1.and.refine(i)) refine(i) = .false.
      end do


! If no leaf blocks were marked for refinement then return.
 21   if(icontinue.eq.0) return


!!!!!!!!!!      
! CHECK FOR neighboring blocks which differ by more than 1 level of refinement
      
      do i = 1,lnblocks
        do iz = 1,1+k3d
          do iy = 1,1+k2d
            do ix = 1,2
              ref_test(ix,iy,iz,i) = .FALSE.
            end do
          end do
        end do
      end do

#ifdef DEBUG
      write(*,*) 'amr_check_refine : proc ',mype,' step 1', & 
     &     ' refine(1:lnblocks) ',refine(1:lnblocks)
#endif /* DEBUG */
! Step 1 - set up an array, ref_test, which contains the 
!          logical refine flag for each child of every local block

      nrecv = 0
      do i = 1,lnblocks
         do j = 1,nchild
            if (child(1,j,i).gt.0) then

              if (j.eq.1) then
                 ix = 1
                 iy = 1
                 iz = 1
              elseif (j.eq.2) then
                 ix = 2
                 iy = 1
                 iz = 1
              elseif (j.eq.3) then
                 ix = 1
                 iy = 2
                 iz = 1
              elseif (j.eq.4) then
                 ix = 2
                 iy = 2
                 iz = 1
              elseif (j.eq.5) then
                 ix = 1
                 iy = 1
                 iz = 2
              elseif (j.eq.6) then
                 ix = 2
                 iy = 1
                 iz = 2
              elseif (j.eq.7) then
                 ix = 1
                 iy = 2
                 iz = 2
              elseif (j.eq.8) then
                 ix = 2
                 iy = 2
                 iz = 2
              endif

               if (child(2,j,i).ne.mype) then
                  nrecv = nrecv + 1
                  call MPI_logical_IRECV(ref_test(ix,iy,iz,i), & 
     &                           1, & 
     &                           MPI_LOGICAL, & 
     &                           child(2,j,i), & 
     &                           child(1,j,i), & 
     &                           amr_mpi_meshComm, & 
     &                           reqr(nrecv), & 
     &                           ierr)
               else
                  ref_test(ix,iy,iz,i) = refine(child(1,j,i))
               end if
            end if
         end do
      end do
           
      nsend = 0
      do i = 1,lnblocks
            ipar = parent(1,i)
            ipar_proc = parent(2,i)
            if (ipar.ge.1) then
               if (ipar_proc.ne.mype) then
                  nsend = nsend + 1
                  call MPI_logical_SSEND(refine(i), & 
     &                           1, & 
     &                           MPI_LOGICAL, & 
     &                           ipar_proc, & 
     &                           i, & 
     &                           amr_mpi_meshComm, & 
     &                           ierr)
               end if
            end if
      end do

      if (nrecv.gt.0) then
        call MPI_WAITALL(nrecv,reqr,statr,ierr)
      end if


! Step 2 - now cycle over all faces, and for each set up the 
!          array, ref_testt, which contains the logical refine flag
!          for each child of the neighbor block across this face.

      do j = 1,nfaces     

#ifdef DEBUG
      Call MPI_BARRIER(amr_mpi_meshComm, ierr)
      write(*,*) 'amr_check_refine : proc ',mype,' step 2',' jface ', & 
     &            j
#endif /* DEBUG */
         do i = 1,lnblocks
            ref_testt(:,:,:,i) = .FALSE.

      if(spherical_pm) then
            ref_testts(:,:,:,i) = .FALSE.
      endif

         end do

         if (j.eq.1) then
            jr = 2
         elseif (j.eq.2) then
            jr = 1
         elseif (j.eq.3) then
            jr = 4
         elseif (j.eq.4) then
            jr = 3
         elseif (j.eq.5) then
            jr = 6
         elseif (j.eq.6) then
            jr = 5
         end if
        
         nrecv = 0
         do i = 1,lnblocks

! do receive from opposite index face, unless opposite face is
! a theta face in spherical coords
         lrecv=.true.

      if(spherical_pm) then
         if(j.eq.3.and.abs(bnd_box(2,2,i)-pi).lt.eps) lrecv=.false.
         if(j.eq.4.and.abs(bnd_box(1,2,i)).lt.eps) lrecv=.false.
      endif

            if (lrecv) then

            if (neigh(1,jr,i).gt.0) then
               if (neigh(2,jr,i).ne.mype) then
                  nrecv = nrecv + 1
#ifdef DEBUG
      write(*,*) 'amr_check_refine : proc ',mype,' blk ',i, & 
     &           ' nrecv ',nrecv, & 
     &           ' source pe ',neigh(2,jr,i), & 
     &           ' posting n recv testt tag ',neigh(1,jr,i),jr
#endif /* DEBUG */
                  call MPI_logical_IRECV(ref_testt(1,1,1,i), & 
     &                           2**ndim, & 
     &                           MPI_LOGICAL, & 
     &                           neigh(2,jr,i), & 
     &                           neigh(1,jr,i), & 
     &                           amr_mpi_meshComm, & 
     &                           reqr(nrecv), & 
     &                           ierr)
               end if
            end if

            end if


      if(spherical_pm) then
! additional receives needed by polar blocks
       lrecv=.true.
       jr0 = jr
       if(j.eq.3.and.abs(bnd_box(1,2,i)).lt.eps) jr0 = 3
       if(j.eq.4.and.abs(bnd_box(1,2,i)).lt.eps) lrecv=.false.
       if(j.eq.4.and.abs(bnd_box(2,2,i)-pi).lt.eps) jr0 = 4
       if(j.eq.3.and.abs(bnd_box(2,2,i)-pi).lt.eps) lrecv=.false.
       if(abs(bnd_box(1,2,i)).lt.eps.and. & 
     &        abs(bnd_box(2,2,i)-pi).lt.eps) then
             lrecv=.true.
        endif
            if (lrecv .and. jr0.eq.j) then
            if (neigh(1,jr0,i).gt.0) then
               if (neigh(2,jr0,i).ne.mype) then
                  nrecv = nrecv + 1
#ifdef DEBUG
      write(*,*) 'amr_check_refine : proc ',mype,' blk ',i, & 
     &           ' nrecv ',nrecv, & 
     &           ' source pe ',neigh(2,jr0,i), & 
     &           ' posting recv testts tag ',neigh(1,jr0,i),jr0
#endif /* DEBUG */

! Note the array ref_testts is larger than is really needed. More careful
! coding may reduce the size.
                  call MPI_logical_IRECV(ref_testts(1,1,1,i), & 
     &                           2**ndim, & 
     &                           MPI_LOGICAL, & 
     &                           neigh(2,jr0,i), & 
     &                           neigh(1,jr0,i), & 
     &                           amr_mpi_meshComm, & 
     &                           reqr(nrecv), & 
     &                           ierr)
               end if
            end if
            end if
         end if  ! end if (spherical_pm

         end do
         
         nsend = 0
         do i = 1,lnblocks

           lsend = .true.

           if(lsend) then
               ineigh = neigh(1,j,i)
               ineigh_proc = neigh(2,j,i)
               if (ineigh.ge.1) then
                  if (ineigh_proc.ne.mype) then
                     nsend = nsend + 1
#ifdef DEBUG
      write(*,*) 'amr_refine_block : proc ',mype,' blk ',i, & 
     &          ' nsend ',nsend, & 
     &          ' dest pe ',ineigh_proc, & 
     &          ' posting send ref_test tag ',i,' neigh() ',neigh(:,j,i)
#endif /* DEBUG */
                     call MPI_logical_SSEND(ref_test(1,1,1,i), & 
     &                              2**ndim, & 
     &                              MPI_LOGICAL, & 
     &                              ineigh_proc, & 
     &                              i, & 
     &                              amr_mpi_meshComm, & 
     &                              ierr)
                  else


!                     ref_testt(:,:,:,ineigh) = ref_test(:,:,:,i)
!start of new test section
         lrecv=.true.
      if(spherical_pm) then
         if(j.eq.3.and.abs(bnd_box(2,2,i)-pi).lt.eps) lrecv=.false.
         if(j.eq.4.and.abs(bnd_box(1,2,i)).lt.eps) lrecv=.false.
      endif
            if (lrecv) then
            if (neigh(1,j,i).gt.0) then
               if (neigh(2,j,i).eq.mype) then
                     ref_testt(:,:,:,neigh(1,j,i)) = ref_test(:,:,:,i)
               end if
            end if
            end if

      if(spherical_pm) then
! additional receives needed by polar blocks
       lrecv=.true.
       jr0 = jr
       if(j.eq.3.and.abs(bnd_box(1,2,i)).lt.eps) jr0 = 3
       if(j.eq.4.and.abs(bnd_box(1,2,i)).lt.eps) lrecv=.false.
       if(j.eq.4.and.abs(bnd_box(2,2,i)-pi).lt.eps) jr0 = 4
       if(j.eq.3.and.abs(bnd_box(2,2,i)-pi).lt.eps) lrecv=.false.
       if(abs(bnd_box(1,2,i)).lt.eps.and. & 
     &        abs(bnd_box(2,2,i)-pi).lt.eps) then
             lrecv=.true.
       endif
            if (lrecv .and. jr0.eq.j) then
            if (neigh(1,jr0,i).gt.0) then
               if (neigh(2,jr0,i).eq.mype) then
! Note the array ref_testts is larger than is really needed. More careful
! coding may reduce the size.
                  ref_testts(:,:,:,i) = ref_test(:,:,:,i)
               end if
            end if
            end if
      endif

!end of new test section


                  end if
               end if

           end if

         end do
#ifdef DEBUG
      write(*,*) 'amr_check_refine : proc ',mype,' waiting jface ',j, & 
     &             ' testt nrecv ',nrecv,' nsend ',nsend
#endif /* DEBUG */
         
         if (nrecv.gt.0) then
            call MPI_WAITALL(nrecv,reqr,statr,ierr)
         end if

#ifdef DEBUG
      Call MPI_BARRIER(amr_mpi_meshComm, ierr)
      write(*,*) 'amr_check_refine : proc ',mype,' step 3',' jface ', & 
     &            j
#endif /* DEBUG */

! Step 3 - if any elements of ref_testt are true, then make sure
!          that the elements of ref_test across face j are also
!          set to true. (Note: j here is the face as indexed by
!          the neighbor of block i.)
         
         do i = 1,lnblocks
           
            if ( any(ref_testt(1:2,:,:,i)) & 
     &         .or.(spherical_pm .and. any(ref_testts(:,:,:,i))) & 
     &         )  then

               
               if (j.eq.1) then
                  
                  ix = 2
                  ix_n = 1
                  do iz = 1,1+k3d
                     do iy = 1,1+k2d
                        ref_test(ix,iy,iz,i) =  & 
     &                       ref_test(ix,iy,iz,i) .or.  & 
     &                       ref_testt(ix_n,iy,iz,i)
                     end do
                  end do
                  
               elseif (j.eq.2) then
                  
                  ix = 1
                  ix_n = 2
                  do iz = 1,1+k3d
                     do iy = 1,1+k2d
                        ref_test(ix,iy,iz,i) =  & 
     &                       ref_test(ix,iy,iz,i) .or.  & 
     &                       ref_testt(ix_n,iy,iz,i)
                     end do
                  end do
              
               elseif (j.eq.3) then
                  
                  iy = 2
                  iy_n = 1          
                  do iz = 1,1+k3d
                     do ix = 1,2
                        ref_test(ix,iy,iz,i) =  & 
     &                       ref_test(ix,iy,iz,i) .or.  & 
     &                       ref_testt(ix,iy_n,iz,i) & 
     &           .or. (spherical_pm .and. ref_testts(ix,iy,iz,i))
                     end do
                  end do
                  
               elseif (j.eq.4) then
                  
                  iy = 1
                  iy_n = 2
                  do iz = 1,1+k3d
                     do ix = 1,2
                        ref_test(ix,iy,iz,i) =  & 
     &                       ref_test(ix,iy,iz,i) .or.  & 
     &                       ref_testt(ix,iy_n,iz,i) & 
     &           .or. (spherical_pm .and. ref_testts(ix,iy,iz,i))
                     end do
                  end do
                  
               elseif (j.eq.5) then
                  
                  iz = 2
                  iz_n = 1
                  do iy = 1,1+k2d
                     do ix = 1,2
                        ref_test(ix,iy,iz,i) =  & 
     &                       ref_test(ix,iy,iz,i) .or.  & 
     &                       ref_testt(ix,iy,iz_n,i)
                     end do
                  end do
                  
               elseif (j.eq.6) then
                  
                  iz = 1
                  iz_n = 2
                  do iy = 1,1+k2d
                     do ix = 1,2
                        ref_test(ix,iy,iz,i) =  & 
     &                       ref_test(ix,iy,iz,i) .or.  & 
     &                       ref_testt(ix,iy,iz_n,i)
                     end do
                  end do
                  
               end if

            end if
            
         end do
      end do

#ifdef DEBUG
      write(*,*) 'amr_check_refine : proc ',mype,' step 4'
#endif /* DEBUG */
! Step 4      
! SET REFINE FLAGS BASED ON ref_test flags and repeat refinement process
! if necessary


      repeat = .FALSE.
      do i = 1,lnblocks

         do iz = 1,1+k3d
           do iy = 1,1+k2d
             do ix = 1,2
           
               if (ref_test(ix,iy,iz,i).and.nodetype(i).eq.1 & 
     &               .and..not.refine(i)) then
                 repeat = .TRUE.
                 refine(i) = .TRUE.
                 derefine(i) = .FALSE.
               else if (ref_test(ix,iy,iz,i)) then
                  ! Setting refine flag here results in cancellation of derefine on
                  ! children of this block when amr_check_derefine processes flags.
                  ! But note that turning off of the refine flags on all non-leaf blocks
                  ! still remains to be done before amr_refine_blocks or amr_derefine_blocks
                  ! get called. - KW
                  refine(i) = .TRUE.
                  derefine(i) = .FALSE.
               end if

             end do
           end do
         end do
               
      end do

! cycle through all processors to see if any repeat, if so all repeat
! this should be done via a scan function

      lreduce_datain(1) = repeat
      call mpi_logical_allreduce (lreduce_datain,lreduce_dataout, & 
     &                    1,MPI_LOGICAL, & 
     &                    MPI_LOR,amr_mpi_meshComm,ierr)
      repeat_t = lreduce_dataout(1)
      repeat = repeat_t

      if (repeat) then
         go to 21
      end if
#ifdef DEBUG_FLOW_TRACE
         write(*,*) 'pe ',mype,' exiting amr_check_refine'
#endif /* DEBUG_FLOW_TRACE */

      return
      end subroutine amr_check_refine
