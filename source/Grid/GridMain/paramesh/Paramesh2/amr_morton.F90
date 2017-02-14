      subroutine amr_morton_order (lnblocks_old,new_loc,nprocs,mype)


! $RCSfile: amr_morton.F90,v $
! $Revision: 1.2 $
! $Date: 2004/09/07 00:58:09 $


! By K. Olson (NASA/GSFC and GMU) 11/96

use physicaldata
      use tree
      implicit none
      include 'mpif.h'


      integer nprocs,mype

      real tot_work
      real time_exe

      integer i,j,minblocks
      integer lnblocks_old,tot_blocks
      integer iii,lb,ierr
      integer mort_no(2*maxblocks_tr)
      integer new_loc(2,maxblocks_tr)
      integer old_loc(2,maxblocks_tr)
      logical :: first = .TRUE.
      save first

      integer, parameter :: FAIL = -1


! compute morton numbers for each cell

      call amr_compute_morton (mort_no)
         
! Sort these morton numbers into order. The subroutine amr_sort_morton
! returns the array new_loc which gives the new locations that 
! the cells are to move to (local address is the first arguement
! and processor number is the second).

      new_loc(:,:) = -1

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file='amr_log',status='unknown', & 
     &        position='append')
         write (30,*) ' calling AMR_SORT_MORTON '
         print *, ' calling AMR_SORT_MORTON '
         close(30)
      end if
#endif

      call amr_sort_morton (mort_no,new_loc,old_loc,nprocs,mype)

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file='amr_log',status='unknown', & 
     &        position='append')
         write (30,*) ' time sort_morton = ',time_exe
         print *,' time sort_morton = ',time_exe
         close(30)
      end if
#endif
      first = .FALSE.
         
! The following call to sort_by_work attempts to realign the 
! sorted list returned by sort_morton such that the work load is 
! balanced across processors.

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file='amr_log',status='unknown', & 
     &        position='append')
         write (30,*) ' calling SORT_BY_WORK '
         print *, ' calling SORT_BY_WORK '
         close(30)
      end if
#endif

      call MPI_ALLREDUCE (lnblocks,tot_blocks,1,MPI_INTEGER, & 
     &                    MPI_SUM,MPI_COMM_WORLD,ierr)

      if (tot_blocks.gt.2*nprocs) then
         call amr_sort_by_work (new_loc,nprocs,mype)
      end if

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file='amr_log',status='unknown', & 
     &        position='append')
         write (30,*) ' time sort_by_work = ',time_exe
         print *,' time sort_by_work = ',time_exe
         close(30)
      end if

      if (mype.eq.0) then
         open (unit=30,file='amr_log',status='unknown', & 
     &        position='append')
         write (30,*) ' calling MIGRATE_TREE '
         print *, ' calling MIGRATE_TREE '
         close(30)
      end if
#endif

      call amr_migrate_tree_data (new_loc,nprocs,mype)

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file='amr_log',status='unknown', & 
     &        position='append')
         write (30,*) ' time migrate_tree = ',time_exe
         print *,' time migrate_tree = ',time_exe
         close(30)
      end if
#endif

! 2) move blocks of data to new locations

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file='amr_log',status='unknown', & 
     &        position='append')
         write (30,*) ' CALLING redist_blk '
         print *,' CALLING redist_blk '
         close(30)
      end if
#endif

      call amr_redist_blk(new_loc,nprocs,mype,lnblocks_old)

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file='amr_log',status='unknown', & 
     &        position='append')
         write (30,*) ' time redist_blks = ',time_exe
         print *,' time redist_blks = ',time_exe
         close(30)
      end if
#endif

      lnblocks = new_lnblocks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file='amr_log',status='unknown', & 
     &        position='append')
         write (30,*) ' DONE MORTON ORDERING'
         write (30,*) ' '
         print *, ' DONE MORTON ORDERING '
         print *,' '
         close(30)
      end if
#endif

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine amr_compute_morton (mort_no)

! This subroutine computes the morton numbers of each cell by interleaving 
! bits in x, y, then z order

! Returns -> mort_no

! By K. Olson (NASA/GSFC and GMU) 12/96

use physicaldata
      use tree
      implicit none
      include 'mpif.h'


      real    xmin,ymin,zmin,xmin_loc,ymin_loc,zmin_loc
      real    size_min_loc,size_min
      integer ix(maxblocks_tr),iy(maxblocks_tr),iz(maxblocks_tr)
      integer mort_no(2*maxblocks_tr),ipos,ipos2,i,j,k
      integer max_level,max_level_loc,ierr

      real xyz_loc_vector(3), xyz_min_loc_vector(3)

! 1) find minimum size of meshes

      size_min_loc = 1.e10
      do i = 1,lnblocks

         if (nodetype(i).eq.1) then

            do j = 1,ndim
               size_min_loc = min(size_min_loc,bsize(j,i))
            end do
            
         end if

      end do

! 1.1) find global size_min across processors
      
      call MPI_ALLREDUCE (size_min_loc,size_min,1, & 
     &                    MPI_DOUBLE_PRECISION, & 
     &                    MPI_MIN,MPI_COMM_WORLD,ierr)

! 2) find local minimum values of x, y, and z

      xmin_loc = 1.e10
      ymin_loc = 1.e10
      zmin_loc = 1.e10
      do i = 1,lnblocks

         if (nodetype(i).eq.1) then

            xmin_loc = min(coord(1,i)-(bsize(1,i)/2.),xmin_loc)
            if (ndim.ge.2) then
               ymin_loc = min(coord(2,i)-(bsize(2,i)/2.),ymin_loc)
            end if
            if (ndim.eq.3) then
               zmin_loc = min(coord(3,i)-(bsize(3,i)/2.),zmin_loc)
            end if

         end if

      end do

! 2.1) find global min^s across processors


! old method -- this is slow
!      call MPI_ALLREDUCE (xmin_loc,xmin,1,MPI_DOUBLE_PRECISION,
!     &                    MPI_MIN,MPI_COMM_WORLD,ierr)

!      if (ndim.ge.2) then
!         call MPI_ALLREDUCE (ymin_loc,ymin,1,MPI_DOUBLE_PRECISION,
!     &                       MPI_MIN,MPI_COMM_WORLD,ierr)
!      end if
!
!      if (ndim.eq.3) then
!         call MPI_ALLREDUCE (zmin_loc,zmin,1,MPI_DOUBLE_PRECISION,
!     &                       MPI_MIN,MPI_COMM_WORLD,ierr)
!      end if


! pack the 3 allreduces into a single allreduce with a vector of the 
! minimum in each coordinate
      xyz_loc_vector(1) = xmin_loc
      xyz_loc_vector(2) = ymin_loc
      xyz_loc_vector(3) = zmin_loc

! reduce in all ndim dimensions
      call MPI_Allreduce(xyz_loc_vector, xyz_min_loc_vector, ndim, & 
     &     MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)

! unpack the minimums
      xmin = xyz_min_loc_vector(1)
      if (ndim .ge. 2) ymin = xyz_min_loc_vector(2)
      if (ndim .ge. 3) zmin = xyz_min_loc_vector(3)
      

! 3) compute ix,iy, and iz

      do i = 1,lnblocks
         
!         ix(i) = int((coord(1,i)-xmin)/size_min)
!         iy(i) = int((coord(2,i)-ymin)/size_min)
!         if (ndim.eq.3) then
!          iz(i) = int((coord(3,i)-zmin)/size_min)
!         else
!          iz(i) = 0
!         end if

!  CODE THAT FOLLOWS gives morton ordering according to a true tree 
!  structure. I.e. child morton numbrs are at least equal to or
!  greater than their parents number.

         ix(i) = int((coord(1,i)-xmin)/bsize(1,i))
         if (ndim.ge.2) then
            iy(i) = int((coord(2,i)-ymin)/bsize(2,i))
         else
            iy(i) = 0
         end if
         if (ndim.eq.3) then
            iz(i) = int((coord(3,i)-zmin)/bsize(3,i))
         else
            iz(i) = 0
         end if
           
      end do
      
! 4) now interleave bits of ix, iy, and iz to get the morton numbers

      do i = 1,maxblocks_tr*2
         mort_no(i) = -1
      enddo
      
      do i = 1,lnblocks
        
 
         ipos = 0
         ipos2 = 0
         mort_no(i) = 0
         
         do while (ipos.lt. 31-3)
              
            call mvbits (ix(i),ipos2,1,mort_no(i),ipos)
            ipos = ipos + 1
            if(ndim.ge.2) then
               call mvbits (iy(i),ipos2,1,mort_no(i),ipos)
               ipos = ipos + 1
            endif
            if(ndim.eq.3) then
               call mvbits (iz(i),ipos2,1,mort_no(i),ipos)
               ipos = ipos + 1
            endif
            
            ipos2 = ipos2 + 1
              
         end do

      end do

!  CODE THAT FOLLOWS gives morton ordering according to a true tree 
!  structure. I.e. child morton numbers are at least equal to or
!  greater than their parents number.
      
! 5.1) find maximum level

      max_level_loc = 0
      do i = 1,lnblocks

            max_level_loc = max(lrefine(i),max_level_loc)
            
      end do

      call MPI_ALLREDUCE (max_level_loc,max_level,1,MPI_INTEGER, & 
     &                    MPI_MAX,MPI_COMM_WORLD,ierr)

! 5) now shift bits to the left by max_levels - level

      do i = 1,lnblocks

!         mort_no(i) = ishft(mort_no(i),3*(max_level-lrefine(i)))
         mort_no(i) = ishft(mort_no(i),ndim*(max_level-lrefine(i)))
         
      end do

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine amr_sort_morton (mort_no,new_loc,old_loc,nprocs,mype)

! Subroutine to sort morton numbers

! Input -> vector of morton numbers, nprocs (no. of processors)

! Output -> new locations that cells are to migrate to
!           new_loc(1,i) is local id to move cell i to
!           new_loc(2,i) is processor id to move cell i to

! Sorting is done without regard to work here.  The new_loc^s returned
! are computed assuming equal (or nearly so) numbers of cells per
! processor.

! By K. Olson (NASA/GSFC and GMU) 11/96

use physicaldata
      use tree
      implicit none
      include 'mpif.h'


      real unkt(nvar,1:iu_bnd,1:ju_bnd,1:ku_bnd)
      integer nprocs,mype
      integer lnblocks2,tot_blocks,no_per_proc,idi,idp
      integer new_loc(2,maxblocks_tr)
      integer old_loc(2,maxblocks_tr)
      integer mort_no(2*maxblocks_tr),i,j,k
      integer irnkg(2*maxblocks_tr)
      integer excess,nprocs_y,nprocs_x,irnkg_s
      integer ierr,next_loc
      integer ix(2*maxblocks_tr),lnblocks_left
      integer is,ie,js,je,ks,ke

!      call amr_bi_sort(mort_no,irnkg,lnblocks)

      do i = 1,lnblocks
         ix(i) = i
      end do
      if (lnblocks.gt.0) then
         call amr_int_heapsort2(mort_no,ix,lnblocks) 
      endif
      do i = 1,lnblocks
         irnkg(ix(i)) = i
      end do

! TEST
      lnblocks_left = 0
      call MPI_SCAN (lnblocks,lnblocks_left,1, & 
     &               MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD, & 
     &               ierr)
      lnblocks_left = lnblocks_left - lnblocks
      do j = 1,lnblocks
        irnkg(j) = irnkg(j) + lnblocks_left
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! 1) Compute total list length.

! I copy lnblocks to lnblocks2 since lnblocks2 can be put in a save statement.
      lnblocks2 = lnblocks 
      call MPI_ALLREDUCE (lnblocks2,tot_blocks,1,MPI_INTEGER, & 
     &                    MPI_SUM,MPI_COMM_WORLD,ierr)

      no_per_proc = tot_blocks/nprocs

      excess = tot_blocks - no_per_proc*nprocs
      nprocs_y = (no_per_proc+1)*nprocs - tot_blocks
! no. of processors which will get no_per_proc + 1 blocks
      nprocs_x = nprocs - nprocs_y
! rank in list which divides those which go on processor with one number
! of blocks from those which go on another set of blocks w. a different
! no. of blocks
      irnkg_s = nprocs_x*(no_per_proc+1)

! 2) Compute new_locs from rankings (irnkg) returned by amr_bi_sort.
!    The following divides blocks evenly among processors without regard to
!    work.

      do i = 1,lnblocks

         idp = (irnkg(i)-1)/(no_per_proc+1) ! processor to send to
         if (irnkg(i).le.irnkg_s) then
            idi = mod((irnkg(i)-1),no_per_proc+1) + 1 ! rank inside 
                                                      ! local array
                                                      ! to write to
               
         else
            idp = (irnkg(i)-irnkg_s-1)/(no_per_proc) ! processor to send to
            idp = idp + nprocs_x
            idi = mod((irnkg(i)-irnkg_s-1),no_per_proc) + 1 ! rank inside 
                                                            ! local array
                                                            ! to write to
         end if

         new_loc(1,i) = idi
         new_loc(2,i) = idp
!         new_loc(1,i) = irnkg(i)
!         new_loc(2,i) = mype
          
      end do

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine amr_sort_by_work (new_loc,nprocs,mype)

! Subroutine to balance work load

! on input takes list sorted by morton ordering
! on output returns new values of new_loc

! By K. Olson (NASA/GSFC and GMU) 11/96

use physicaldata
      use tree
      implicit none
      include 'mpif.h'



      integer nprocs,mype,ierr,errorcode
      
      integer lnblocks2,lnblocksl

      real temp

      real work(maxblocks_tr),workt(maxblocks_tr),loc_work,tot_work
      real work_per_proc,work_left ! dim. set to give max. no. of
                                         ! procs expected

      integer i,j,k,idi
      integer itemp,itemp2
      integer pid(maxblocks_tr),lid(maxblocks_tr),lid2(maxblocks_tr)
      integer pidt,lidt,lid_old,tot_blocks
      integer left,right
      integer nsend,nrecv,max_nsend,max_nrecv
      integer statr(MPI_STATUS_SIZE,2*maxblocks_tr)
      integer stat (MPI_STATUS_SIZE)
      integer reqr(2*maxblocks_tr)
      
      integer new_loc(2,maxblocks_tr)
      integer new_loc_temp(2,maxblocks_tr)
      integer old_loc(2,maxblocks_tr)

      logical repeat,repeatt

! initialize work arrary and temp work arrary

      do i = 1,maxblocks_tr

         work(i) = -1.
         workt(i) = -1.

      end do

! assign values to work array

      do i = 1,maxblocks_tr

         work(i) = 0.
         work(i) = work_block(i)
!         if (nodetype(i).eq.1) work(i) = 2.
!         if (nodetype(i).eq.2) work(i) = 1.

      end do

! move work to temp work array workt

      call fill_old_loc (new_loc,old_loc,nprocs,mype)

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30, & 
     &        file='amr_log', & 
     &        position='append', & 
     &        status='unknown', & 
     &        form='formatted')
         write (30,*) ' DONE fill_old_loc '
         close(30)
      end if
#endif

      lnblocks2 = 0
      do i = 1,maxblocks_tr
        if (old_loc(1,i).gt.-1) then
          lnblocks2 = lnblocks2 + 1
        end if
      end do

      nrecv = 0
      do i = 1,lnblocks2
        if (old_loc(2,i).ne.mype) then
          nrecv = nrecv + 1
          call MPI_IRECV(workt(i),1,MPI_DOUBLE_PRECISION, & 
     &         old_loc(2,i),old_loc(1,i),MPI_COMM_WORLD, & 
     &         reqr(nrecv),ierr)
        end if
      end do

      nsend = 0
      do i = 1,lnblocks

        if (new_loc(2,i).ne.mype) then
          nsend = nsend + 1
          call MPI_SSEND(work(i),1,MPI_DOUBLE_PRECISION, & 
     &         new_loc(2,i),i,MPI_COMM_WORLD,ierr)
        else
          workt(new_loc(1,i)) = work(i)
        end if

      end do
         
      if (nrecv.gt.0) then
        call MPI_WAITALL(nrecv,reqr,statr,ierr)
      end if



      do i = 1,lnblocks2
         work(i) = workt(i)
      end do

! SUM total work within each processosr

      if (lnblocks2.gt.0) then
        workt(1) = work(1)
      else
        workt(1) = 0
      end if
      do i = 2,lnblocks2

         workt(i) = workt(i-1) + work(i)

      end do

! SUM work across processors

      if (lnblocks2.gt.0) then
        loc_work = workt(lnblocks2)
      else
        loc_work = 0
      end if

      call MPI_ALLREDUCE (loc_work,tot_work,1,MPI_DOUBLE_PRECISION, & 
     &                    MPI_SUM,MPI_COMM_WORLD,ierr)

! Compute work per processor

      work_per_proc = tot_work/nprocs


! Compute final work by looking left

      work_left = 0.

      call MPI_SCAN (loc_work,work_left,1, & 
     &               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, & 
     &               ierr)

      work_left = work_left - loc_work
      do j = 1,lnblocks2
         workt(j) = workt(j) + work_left
      end do

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30, & 
     &        file='amr_log', & 
     &        position='append', & 
     &        status='unknown', & 
     &        form='formatted')
         write (30,*) ' DONE compute work_left '
         close(30)
      end if
#endif
! compute processor ids

      do i = 1,maxblocks_tr

         pid(i) = 0
         lid(i) = 0

      end do

      do i = 1,lnblocks2

         pid(i) = int((workt(i)-1.)/work_per_proc)
         if (pid(i).lt.0) pid(i) = 0
         if (pid(i).gt.nprocs-1) pid(i) = nprocs-1

      end do

! compute local ids
      
      lid(1) = 1
      do i = 2,lnblocks2

         lid(i) = lid(i-1) + 1
         if (pid(i-1).lt.pid(i)) lid(i) = 1  ! start a new group

      end do

      do i = 1,maxblocks_tr
         lid2(i) = lid(i)
      end do

      left = mype - 1
      right = mype + 1
      if (mype.eq.0) left = MPI_PROC_NULL
      if (mype.eq.nprocs-1) right = MPI_PROC_NULL

      pidt = 0

      call MPI_SENDRECV  & 
     &     (pid(lnblocks2),1,MPI_INTEGER,right,1, & 
     &      pidt,          1,MPI_INTEGER,left, 1, & 
     &      MPI_COMM_WORLD,stat,ierr)

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30, & 
     &        file='amr_log', & 
     &        position='append', & 
     &        status='unknown', & 
     &        form='formatted')
         write (30,*) ' STARTING loop 27 '
         close(30)
      end if
#endif

      lidt = 0
 27   lid_old = lidt ! lid_old stores last fetched value of lid to left

      lidt = 0

      call MPI_SENDRECV  & 
     &     (lid(lnblocks2),1,MPI_INTEGER,right,1, & 
     &      lidt,          1,MPI_INTEGER,left, 1, & 
     &      MPI_COMM_WORLD,stat,ierr)

      do j = 1,lnblocks2
         
         if (pidt.eq.pid(j)) then ! if pidt (which was fetched)
                                  ! equals local pid then the list
                                  ! has been split across processors
               
            lid(j) = lid2(j) + lidt
            
         end if

      end do
      
      repeat = .FALSE.
      if (lidt.ne.lid_old) repeat = .TRUE.
      call MPI_ALLREDUCE (repeat,repeatt,1,MPI_LOGICAL, & 
     &                    MPI_LOR,MPI_COMM_WORLD,ierr)
      if (repeatt) go to 27

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30, & 
     &        file='amr_log', & 
     &        position='append', & 
     &        status='unknown', & 
     &        form='formatted')
         write (30,*) ' DONE loop 27 '
         close(30)
      end if
#endif

! now reorder according to new pid and lid numbers

      nrecv = 0
      do i = 1,lnblocks
         if (new_loc(2,i).ne.mype) then
            nrecv = nrecv + 1
            call MPI_IRECV(new_loc_temp(1,i),1,MPI_INTEGER, & 
     &           new_loc(2,i),new_loc(1,i),MPI_COMM_WORLD, & 
     &           reqr(nrecv),ierr)
         else
            new_loc_temp(1,i) = lid(new_loc(1,i))
         end if
      end do

      nsend = 0
      do i = 1,lnblocks2
        if (old_loc(2,i).ne.mype) then
           nsend = nsend + 1
           call MPI_SSEND(lid(i),1,MPI_INTEGER, & 
     &          old_loc(2,i),i, & 
     &          MPI_COMM_WORLD,ierr)
        end if
      end do

      if (nrecv.gt.0) then
        call MPI_WAITALL(nrecv,reqr,statr,ierr)
      end if

      nrecv = 0
      do i = 1,lnblocks
         if (new_loc(2,i).ne.mype) then
            nrecv = nrecv + 1
            call MPI_IRECV(new_loc_temp(2,i),1,MPI_INTEGER, & 
     &           new_loc(2,i),new_loc(1,i), & 
     &           MPI_COMM_WORLD, & 
     &           reqr(nrecv),ierr)
         else
            new_loc_temp(2,i) = pid(new_loc(1,i))
         end if
      end do
      
      nsend = 0
      do i = 1,lnblocks2
         if (old_loc(2,i).ne.mype) then
           nsend = nsend + 1
           call MPI_SSEND(pid(i),1,MPI_INTEGER, & 
     &          old_loc(2,i),i, & 
     &          MPI_COMM_WORLD,ierr)
        end if
      end do
         
      if (nrecv.gt.0) then
        call MPI_WAITALL(nrecv,reqr,statr,ierr)
      end if

      new_loc(:,1:lnblocks) = new_loc_temp(:,1:lnblocks)

      do i = 1,lnblocks
         if( new_loc(1,i) > maxblocks ) then
            open (unit=30, & 
     &           file='amr_log', & 
     &           position='append', & 
     &           status='unknown', & 
     &           form='formatted')
            write(30,*)
            write(30,*) '[AMR_SORT_BY_WORK] ERROR: block location exceeds MAXBLOCKS limit'
            write(30,*)
            write(30,*) 'mype, i ', mype, i
            write(30,*) 'new_loc(1,i) = local address = ',new_loc(1,i)
            write(30,*) 'new_loc(2,i) = processor     = ',new_loc(2,i)
            write(30,*)
            write(30,*) 'Suggestion: increase MAXBLOCKS, number of processors or modify refinement criteria'
            close(30)
            call Driver_abortFlash("[AMR_SORT_BY_WORK] ERROR: block location exceeds MAXBLOCKS limit")
         endif         
         
         if ( (new_loc(2,i) < 0) .OR. (new_loc(2,i) > nprocs) ) then
            open (unit=30, & 
     &           file='amr_log', & 
     &           position='append', & 
     &           status='unknown', & 
     &           form='formatted')
            write(30,*)
            write(30,*) '[AMR_SORT_BY_WORK] ERROR: target processor number out of bounds'
            write(30,*)
            write(30,*) 'mype, i ', mype, i
            write(30,*) 'new_loc(2,i) = processor = ',new_loc(2,i)
            write(30,*)
            write(30,*) 'Suggestion: increase MAXBLOCKS, number of processors or modify refinement criteria'
            close(30)
            call Driver_abortFlash("[AMR_SORT_BY_WORK] ERROR: target processor number out of bounds")
         endif         

      end do

      new_loc(:,lnblocks+1:maxblocks_tr) = -1

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30, & 
     &        file='amr_log', & 
     &        position='append', & 
     &        status='unknown', & 
     &        form='formatted')
         write (30,*) ' DONE sort_by_work '
         close(30)
      end if
#endif

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine amr_migrate_tree_data (new_loc,nprocs,mype)

! Subroutine move tree data and reconnect all pointers given new_loc

! By K. Olson (NASA/GSFC and GMU) 11/96

use physicaldata
      use tree
      implicit none
      include 'mpif.h'


      integer buf_size
      integer ibuf_size
      parameter(buf_size = mdim+mdim+2*mdim)
      real buffer(buf_size)
      real buffert(buf_size,maxblocks_tr)

      parameter(ibuf_size = 2*mfaces+2*mchild+2+2)
      integer ibuffer(ibuf_size)
      integer ibuffert(ibuf_size,maxblocks_tr)

      logical newchildt(maxblocks_tr)

      integer nprocs,mype
      
      integer neight(2,mfaces,maxblocks_tr)
      integer childt(2,mchild,maxblocks_tr)
      integer parentt(2,maxblocks_tr),lrefinet(maxblocks_tr)
      integer nodetypet(maxblocks_tr)
      integer emptyt(maxblocks_tr)

      real coordt(mdim,maxblocks_tr),sizet(mdim,maxblocks_tr)
      real bnd_boxt(2,mdim,maxblocks_tr)

      real temp,tempb(2)
      real rnum

      integer i,j,k,jj,idi,idp,lb
      integer itemp,cempty,remote_block,remote_pe
      
      integer new_loc(2,maxblocks_tr)
      integer old_loc(2,maxblocks_tr)

      integer statr(MPI_STATUS_SIZE,maxblocks_tr)
      integer reqr(maxblocks_tr)
      integer ierr,nsend,nrecv,new_blocks_old

      logical ltemp

      integer tot_blocks,new_tot_blocks

      call fill_old_loc (new_loc,old_loc,nprocs,mype)

! count no. of new blocks

      new_lnblocks = 0
      do i = 1,maxblocks_tr
        if (old_loc(1,i).gt.-1) then
          new_lnblocks = new_lnblocks + 1
        end if
      end do

! update pointers to parents, children and neighbors

      parentt(:,1:lnblocks) = parent(:,1:lnblocks)
      childt(:,:,1:lnblocks) = child(:,:,1:lnblocks)
      neight(:,:,1:lnblocks) = neigh(:,:,1:lnblocks)

      nrecv = 0
      do i = 1,lnblocks
         if (parent(1,i).gt.0) then
           if (parent(2,i).ne.mype) then
             nrecv = nrecv + 1
             call MPI_IRECV(parentt(1,i),2,MPI_INTEGER, & 
     &            parent(2,i),i,MPI_COMM_WORLD, & 
     &            reqr(nrecv),ierr)
           else
             parentt(:,i) = new_loc(:,parent(1,i))
           end if
         end if
       end do
       
       nsend = 0
       do i = 1,lnblocks
         do j = 1,nchild
           if (child(1,j,i).gt.0) then
             if (child(2,j,i).ne.mype) then
               ! parent is sending to all its children
               nsend = nsend + 1
               call MPI_SSEND (new_loc(1,i),2,MPI_INTEGER, & 
     &              child(2,j,i),child(1,j,i),MPI_COMM_WORLD, & 
     &              ierr)
             end if
           end if
         end do
       end do

      if (nrecv.gt.0) then
        call MPI_WAITALL(nrecv,reqr,statr,ierr)
      end if

      nrecv = 0
      do i = 1,lnblocks
        do j = 1,nchild
          if (child(1,j,i).gt.0) then
            if (child(2,j,i).ne.mype) then
              nrecv = nrecv + 1
              call MPI_IRECV(childt(1,j,i),2,MPI_INTEGER, & 
     &             child(2,j,i),child(1,j,i),MPI_COMM_WORLD, & 
     &             reqr(nrecv),ierr)
            else
              childt(:,j,i) = new_loc(:,child(1,j,i))
            end if
          end if
        end do
       end do
       
       nsend = 0
       do i = 1,lnblocks
         if (parent(1,i).gt.0) then
           if (parent(2,i).ne.mype) then
! child is sending to its parent
             nsend = nsend + 1
             call MPI_SSEND (new_loc(1,i),2,MPI_INTEGER, & 
     &            parent(2,i),i,MPI_COMM_WORLD, & 
     &            ierr)
           end if
         end if
       end do

      if (nrecv.gt.0) then
        call MPI_WAITALL(nrecv,reqr,statr,ierr)
      end if

      nrecv = 0
      do i = 1,lnblocks
        do j = 1,nfaces
          if (neigh(1,j,i).gt.0) then
            if (neigh(2,j,i).ne.mype) then
              nrecv = nrecv + 1
              call MPI_IRECV(neight(1,j,i),2,MPI_INTEGER, & 
     &             neigh(2,j,i),neigh(1,j,i),MPI_COMM_WORLD, & 
     &             reqr(nrecv),ierr)
            else
              neight(:,j,i) = new_loc(:,neigh(1,j,i))
            end if
          end if
        end do
      end do
      
      nsend = 0
      do i = 1,lnblocks
        do j = 1,nfaces
          if (neigh(1,j,i).gt.0) then
            if (neigh(2,j,i).ne.mype) then
              nsend = nsend + 1
              call MPI_SSEND (new_loc(1,i),2,MPI_INTEGER, & 
     &             neigh(2,j,i),i,MPI_COMM_WORLD, & 
     &             ierr)
            end if
          end if
        end do
      end do

      if (nrecv.gt.0) then
        call MPI_WAITALL(nrecv,reqr,statr,ierr)
      end if

      parent(:,1:lnblocks) = parentt(:,1:lnblocks)
      child(:,:,1:lnblocks) = childt(:,:,1:lnblocks)
      neigh(:,:,1:lnblocks) = neight(:,:,1:lnblocks)

! initialize temp buffer array

      do i = 1,maxblocks_tr
         buffert(:,i) = -1.
         ibuffert(:,i) = -1
         newchildt(i) = .FALSE.
      end do

      nrecv = 0
      do i = 1,new_lnblocks
        if (old_loc(2,i).ne.mype) then
          nrecv = nrecv + 1
          call MPI_IRECV(buffert(1,i),buf_size,MPI_DOUBLE_PRECISION, & 
     &         old_loc(2,i),i,MPI_COMM_WORLD, & 
     &         reqr(nrecv),ierr)
          nrecv = nrecv + 1
          call MPI_IRECV(ibuffert(1,i),ibuf_size,MPI_INTEGER, & 
     &         old_loc(2,i),i+maxblocks_tr,MPI_COMM_WORLD, & 
     &         reqr(nrecv),ierr)
          nrecv = nrecv + 1
          call MPI_IRECV(newchildt(i),1,MPI_LOGICAL, & 
     &         old_loc(2,i),i+2*maxblocks_tr,MPI_COMM_WORLD, & 
     &         reqr(nrecv),ierr)
        end if
      end do

      nsend = 0
      do i = 1,lnblocks

! pack buffer for sending

        k = 0
        do j = 1,mdim
          k = k + 1
          buffer(k) = coord(j,i)
        end do
        do j = 1,mdim
          do jj = 1,2
            k = k + 1
            buffer(k) = bnd_box(jj,j,i)
          end do
        end do
        do j = 1,mdim
          k = k + 1
          buffer(k) = bsize(j,i)
        end do

        k = 0
        do j = 1,mchild
          do jj = 1,2
            k = k + 1
            ibuffer(k) = child(jj,j,i)
          end do
        end do
        do j = 1,mfaces
          do jj = 1,2
            k = k + 1
            ibuffer(k) = neigh(jj,j,i)
          end do
        end do
        do j = 1,2
          k = k + 1
          ibuffer(k) = parent(j,i)
        end do
        k = k + 1
        ibuffer(k) = lrefine(i)
        k = k + 1
        ibuffer(k) = nodetype(i)

        if (new_loc(2,i).ne.mype) then
          nsend = nsend + 1
          call MPI_SSEND(buffer(1),buf_size,MPI_DOUBLE_PRECISION, & 
     &         new_loc(2,i),new_loc(1,i), & 
     &         MPI_COMM_WORLD,ierr)
          nsend = nsend + 1
          call MPI_SSEND(ibuffer(1),ibuf_size,MPI_INTEGER, & 
     &         new_loc(2,i),new_loc(1,i)+maxblocks_tr, & 
     &         MPI_COMM_WORLD,ierr)
          nsend = nsend + 1
          call MPI_SSEND(newchild(i),1,MPI_LOGICAL, & 
     &         new_loc(2,i),new_loc(1,i)+2*maxblocks_tr, & 
     &         MPI_COMM_WORLD,ierr)
        else
          buffert(1:buf_size,new_loc(1,i)) = buffer(1:buf_size)
          ibuffert(1:ibuf_size,new_loc(1,i)) = ibuffer(1:ibuf_size)
          newchildt(new_loc(1,i)) = newchild(i)
        end if

      end do

      if (nrecv.gt.0) then
        call MPI_WAITALL(nrecv,reqr,statr,ierr)
      end if
  
 ! no unpack the buffer

      do i = 1,maxblocks_tr

        k = 0
        do j = 1,mdim
          k = k + 1
          coord(j,i) = buffert(k,i)
        end do
        do j = 1,mdim
          do jj = 1,2
            k = k + 1
            bnd_box(jj,j,i) = buffert(k,i)
          end do
        end do
        do j = 1,mdim
          k = k + 1
          bsize(j,i) = buffert(k,i)
        end do

        k = 0
        do j = 1,mchild
          do jj = 1,2
            k = k + 1
            child(jj,j,i) = ibuffert(k,i)
          end do
        end do
        do j = 1,mfaces
          do jj = 1,2
            k = k + 1
            neigh(jj,j,i) = ibuffert(k,i)
          end do
        end do
        do j = 1,2
          k = k + 1
          parent(j,i) = ibuffert(k,i)
        end do
        k = k + 1
        lrefine(i) = ibuffert(k,i)
        k = k + 1
        nodetype(i) = ibuffert(k,i)
        k = k + 1
        empty(i) = 0 ! THIS WAS IT !!!!

        newchild(i) = newchildt(i)

      end do
      
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine fill_old_loc (new_loc,old_loc,nprocs,mype)

use physicaldata
      use tree
      implicit none
      include 'mpif.h'


      integer new_loc(2,maxblocks_tr)
      integer old_loc(2,maxblocks_tr)
      integer statr(MPI_STATUS_SIZE,maxblocks_tr)
      integer reqr(maxblocks_tr)
      integer kk(maxblocks_tr),nrecv,nsend, nrecv_old
      integer nsend_to_proc(0:16384), nrecv_pack(0:16384)
      integer ierr,nprocs,mype
      integer i,j,k


! fill `old_loc' (pointer from new block location back to
! its old, unsorted location)

! count no. of receives to post

! 1) count no. of sends on each proc. to all other procs

      nsend_to_proc(:) = 0
      do i = 1,maxblocks_tr
        if (new_loc(1,i).gt.0) then
        if (new_loc(2,i).ne.mype) then ! its a send
          nsend_to_proc(new_loc(2,i)) =  & 
     &    nsend_to_proc(new_loc(2,i)) + 1
        end if
        end if
      end do

! 2) collect data for `this' proc from other procs so that
!    the total no. of receives to post can be computed


      nrecv = 0

! old way
!      do i = 0,nprocs-1
!        call MPI_REDUCE ( nsend_to_proc(i), nrecv, 1, 
!     &                    MPI_INTEGER, MPI_SUM, i, MPI_COMM_WORLD,
!     &                    ierr)
!      end do

! new way -- this scales a lot better
      call MPI_AllReduce(nsend_to_proc, nrecv_pack, nprocs, & 
     &     MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      
      nrecv = nrecv_pack(MyPE)

      old_loc(:,:) = -1
      do i = 1,nrecv
        call MPI_IRECV(kk(i), & 
     &       1, & 
     &       MPI_INTEGER, & 
     &       MPI_ANY_SOURCE, & 
     &       MPI_ANY_TAG, & 
     &       MPI_COMM_WORLD, & 
     &       reqr(i), & 
     &       ierr)
      end do

      nsend = 0
      do i = 1,maxblocks_tr
        if (new_loc(1,i).gt.0) then
        if (new_loc(2,i).ne.mype) then
          nsend = nsend + 1
          call MPI_SSEND(i, & 
     &         1, & 
     &         MPI_INTEGER, & 
     &         new_loc(2,i),    &  ! PE TO SEND TO
     &         new_loc(1,i),    &  ! THIS IS THE TAG
     &         MPI_COMM_WORLD, & 
     &         ierr)
        else
          old_loc(1,new_loc(1,i)) = i
          old_loc(2,new_loc(1,i)) = mype
        end if
        end if
      end do

      if (nrecv.gt.0) then
         call MPI_WAITALL (nrecv, reqr, statr, ierr)
         do i = 1,nrecv
           old_loc(1,statr(MPI_TAG,i)) = kk(i)
           old_loc(2,statr(MPI_TAG,i)) = statr(MPI_SOURCE,i)
         end do
      end if

      call MPI_BARRIER (MPI_COMM_WORLD,ierr) ! NEEDED ????

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

