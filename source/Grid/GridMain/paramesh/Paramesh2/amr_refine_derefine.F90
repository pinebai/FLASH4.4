subroutine amr_refine_derefine

! $RCSfile: amr_refine_derefine.F90,v $
! $Revision: 1.2 $
! $Date: 2004/09/07 00:58:09 $

  use physicaldata
  use Driver_interface, ONLY : Driver_abortFlash
  use tree
  implicit none
  include 'mpif.h'
  


  real times,time_exe


#ifdef TIMINGS
#include "timer.fh"
  integer itimer1,itimer2
#endif

  character (len=80) :: tblk_string
  
  integer nprocs,mype
  integer i,j,k,l
  integer lnblocks_old
  integer ierr,errorcode
  
  logical refinet(maxblocks_tr)
  logical lparref
  
  integer lnblocks2,tot_blocks,icontinue
  integer icontinue_ref,icontinue_deref
  integer icontinue2
  integer min_blocks, max_blocks, tot_blocksa
  integer min_blocks_leaf, max_blocks_leaf, tot_blocksa_leaf, lnblocks_leaf
  integer tot_ref,tot_deref,lb
  integer parent_blk,parent_pe,cempty,nrefs
  integer new_loc(2,maxblocks_tr)
  integer nrefine
  
  
  
  ! string containers to dump name/value pairs 
  ! of tot_blocks, min_blocks, max_blocks into logfile
  
  character(len=32), dimension(3,2) :: block_buff
  character(len=32)                 :: int_to_str
  
  
  nrefine = 0
  
  call MPI_COMM_SIZE (MPI_COMM_WORLD,nprocs,ierr)
  call MPI_COMM_RANK (MPI_COMM_WORLD,mype,ierr)
  
  refinet(1:lnblocks) = refine(1:lnblocks)
  newchild(:) = .FALSE.


! CHECK derefinements and refinements

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! test to see if any refinements have been requested.
  icontinue=0
  icontinue_deref = 0
  if(lnblocks.gt.0) then
     do l = 1,lnblocks
        if(nodetype(l).eq.1.and.refine(l)) then
           icontinue=1
           goto 10
        endif
     enddo
  endif
10 continue
  call MPI_ALLREDUCE (icontinue,icontinue2,1,MPI_INTEGER, & 
       &                      MPI_MAX,MPI_COMM_WORLD,ierr)
  icontinue_ref = icontinue2
  if(icontinue_ref.eq.0) go to 20
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
#ifdef DEBUG_AMR
  if (mype.eq.0) then
     open (unit=30,file='amr_log',status='unknown', & 
          &        position='append')
     write (30,*) ' starting CHECK_REFINE '
     close(30)
     print *, ' starting CHECK_REFINE '
  end if
#endif
  
  call amr_check_refine (nprocs,mype)
  
20 continue
  
#ifdef DEBUG_AMR
  if (mype.eq.0) then
     open (unit=30,file='amr_log',status='unknown', & 
          &        position='append')
     write (30,*) ' time check_refine = ',time_exe
     print *,' time check_refine = ',time_exe
     write (30,*)' starting CHECK_DEREFINE '
     close(30)
     print *,' starting CHECK_DEREFINE '
  end if
#endif
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call amr_check_derefine (mype)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! test to see if any derefinements have been requested
  icontinue=0
  icontinue_deref = 0
  if(lnblocks.gt.0) then
     do l=1,lnblocks
        if(nodetype(l).eq.1.and.derefine(l)) then
           icontinue=1
           goto 101
        endif
     enddo
  endif
101 continue
  call MPI_ALLREDUCE (icontinue,icontinue2,1,MPI_INTEGER, & 
       &     MPI_MAX,MPI_COMM_WORLD,ierr)
  icontinue_deref = icontinue2
  if(icontinue_ref.eq.0.and.icontinue_deref.eq.0) return

      grid_changed = 1

!!$  if (mype.eq.0) then
!!$     call stamp_logfile ("refinement in progress..." )
!!$  end if
        


        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! NOW Actually refine and derefine the mesh

#ifdef DEBUG_AMR
  if (mype.eq.0) then         
     open (unit=30,file='amr_log',status='unknown', & 
          &        position='append')
     print *,' time check_derefine = ',time_exe
     write (30,*)' time check_derefine = ',time_exe
     write (30,*)' starting REFINE_BLOCKS'
     close(30)
     print *,' starting REFINE_BLOCKS'
  end if
#endif
  
  lnblocks_old = lnblocks
  call amr_refine_blocks (new_loc,nprocs,mype)
  
#ifdef DEBUG_AMR
  if (mype.eq.0) then
     open (unit=30,file='amr_log',status='unknown', & 
          &        position='append')
     print *,' time refine_blocks = ',time_exe
     write (30,*)' time refine_blocks = ',time_exe
     write (30,*) ' starting DEREFINE_BLOCKS'
     close(30)
     print *, ' starting DEREFINE_BLOCKS'
  end if
#endif     
  
  call amr_derefine_blocks(lnblocks_old,new_loc,mype)
  
#ifdef DEBUG_AMR
  if (mype.eq.0) then
     open (unit=30,file='amr_log',status='unknown', & 
          &        position='append')
     write (30,*) ' time derefine_blocks = ',time_exe
     print *,' time derefine_blocks = ',time_exe
     close(30)
  end if
#endif
  
  call MPI_ALLREDUCE (lnblocks,tot_blocks,1,MPI_INTEGER, & 
       &                    MPI_SUM,MPI_COMM_WORLD,ierr)
#ifdef DEBUG_AMR
  if (mype.eq.0) then
     open (unit=30,file='amr_log',status='unknown', & 
          &        position='append')
     write(30,*) ' tot_blocks before ',tot_blocks
     close(30)
     print *,' tot_blocks before ',tot_blocks
  end if
  
  ! I copy lnblocks to lnblocks2 since lnblocks2 can be put in a save statement.
  lnblocks2 = lnblocks 
  call MPI_ALLREDUCE (lnblocks2,max_blocks,1,MPI_INTEGER, & 
       &                    MPI_MAX,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE (lnblocks2,min_blocks,1,MPI_INTEGER, & 
       &                    MPI_MIN,MPI_COMM_WORLD,ierr)
  
  if (mype.eq.0) then
     open (unit=30,file='amr_log',status='unknown', & 
          &        position='append')
     write (30,*) ' max_blocks 1',max_blocks
     write (30,*) ' min_blocks 1',min_blocks
     print *, ' max_blocks 1',max_blocks
     print *, ' min_blocks 1',min_blocks
     close(30)
  end if
#endif
  
  ! set work values
  
#ifdef DEBUG_AMR
  if (mype.eq.0) then
     open (unit=30,file='amr_log',status='unknown', & 
          &        position='append')
     write (30,*) ' starting MORTON ORDERING'
     print *, ' starting MORTON ORDERING '
     close(30)
  end if
#endif
  
  work_block(:) = 0.
  do i = 1,lnblocks
     if (nodetype(i).eq.1) work_block(i) = 2.
     if (nodetype(i).eq.2) work_block(i) = 1.
  end do
  call amr_morton_order (lnblocks_old,new_loc,nprocs,mype)
  
  lnblocks_leaf = 0
  
  do i = 1,lnblocks
     if (nodetype(i) == 1) lnblocks_leaf = lnblocks_leaf + 1
  end do
  
  ! I copy lnblocks to lnblocks2 since lnblocks2 can be put in a save statement.
  lnblocks2 = lnblocks 
  call MPI_ALLREDUCE (lnblocks2,tot_blocksa,1,MPI_INTEGER, & 
       &                    MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE (lnblocks2,max_blocks,1,MPI_INTEGER, & 
       &                    MPI_MAX,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE (lnblocks2,min_blocks,1,MPI_INTEGER, & 
       &                    MPI_MIN,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE (lnblocks_leaf,tot_blocksa_leaf,1,MPI_INTEGER, & 
       &                    MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE (lnblocks_leaf,min_blocks_leaf,1,MPI_INTEGER, & 
       &                    MPI_MIN,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE (lnblocks_leaf,max_blocks_leaf,1,MPI_INTEGER, & 
       &                    MPI_MAX,MPI_COMM_WORLD,ierr)
  
  if (mype.eq.0) then
     
     open (unit=30,file='amr_log',status='unknown', & 
          &        position='append')
     write(30,'(3(a,i10))')                                     &
          ' min_blocks ',min_blocks,         &
          ' max_blocks ',max_blocks,         &
          ' tot_blocks ',tot_blocksa
     close (30)
     
     write(*,'(3(a,i10))')                                      &
          ' min_blocks ',min_blocks,         &
          ' max_blocks ',max_blocks,         &
          ' tot_blocks ',tot_blocksa
     
     ! write to string array and pass it to logfile;
     ! array dimensions are hard-coded in declaration
     
     write (block_buff(1, 1), "(A)") "min_blocks"
     write (int_to_str, "(I32)") min_blocks
     write (block_buff(1, 2), "(A)") trim(adjustl(int_to_str))
     
     write (block_buff(2, 1), "(A)") "max_blocks"
     write (int_to_str, "(I32)") max_blocks
     write (block_buff(2, 2), "(A)") trim(adjustl(int_to_str))
     
     write (block_buff(3, 1), "(A)") "tot_blocks"
     write (int_to_str, "(I32)") tot_blocksa
     write (block_buff(3, 2), "(A)") trim(adjustl(int_to_str))
     
!!$     call stamp_logfile(block_buff, 3, 2, "refined all ")

     write (block_buff(1, 1), "(A)") "min_blocks"
     write (int_to_str, "(I32)") min_blocks_leaf
     write (block_buff(1, 2), "(A)") trim(adjustl(int_to_str))
     
     write (block_buff(2, 1), "(A)") "max_blocks"
     write (int_to_str, "(I32)") max_blocks_leaf
     write (block_buff(2, 2), "(A)") trim(adjustl(int_to_str))
     
     write (block_buff(3, 1), "(A)") "tot_blocks"
     write (int_to_str, "(I32)") tot_blocksa_leaf
     write (block_buff(3, 2), "(A)") trim(adjustl(int_to_str))
     
!!$     call stamp_logfile(block_buff, 3, 2, "refined leaf")

  end if

#ifdef DEBUG_AMR
  if (tot_blocksa.ne.tot_blocks) then
     print *,'[AMR_REFINE_DEREFINE] ERROR: tot_blocksa.ne.tot_blocks ', & 
          &        tot_blocksa,tot_blocks
     call Driver_abortFlash("[AMR_REFINE_DEREFINE] ERROR: tot_blocksa.ne.tot_blocks")
  end if
#endif

! Get the nodetypes of all neighbors and children for each block

  if (mype.eq.0)  then
     open (unit=30,file='amr_log',status='unknown', & 
          &        position='append')
     write (30,*) 'CALLING get_tree'
     close(30)
     print *, 'CALLING get_tree'
  end if
  
  
  call get_tree_nodetypes(nprocs,mype)
  
#ifdef DEBUG_AMR
  if (mype.eq.0) then
     open (unit=30,file='amr_log',status='unknown', & 
          &        position='append')
     write (30,*) ' done REFINE_DEREFINE '
     write (30,*) ' '
     print *, ' done REFINE_DEREFINE '
     print *,' '
     close(30)
  end if
#endif
  
  return
end subroutine amr_refine_derefine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_tree_nodetypes (nprocs,mype)
  use Grid_data, ONLY : gr_msgbuffer
  use physicaldata
  use tree
  include 'mpif.h'
  
  
  integer nprocs,mype
  integer i,j,k,isg,neighs,neighr
  integer jsend
  integer ierr
  integer reqr(maxblocks_tr)
  integer statr(MPI_STATUS_SIZE,maxblocks_tr)
  integer stags(maxblocks_tr),svals(maxblocks_tr)
  integer sprocs(maxblocks_tr)
  integer rtags(maxblocks_tr),rvals(maxblocks_tr)
  integer rprocs(maxblocks_tr)
  

      
!!      neigh_type(:,:) = -1
!!      child_type(:,:) = -1

  neighr = 0
  do isg = 1,lnblocks
     do j = 1,nfaces
        if(neigh(1,j,isg).gt.-1) then
           if(neigh(2,j,isg).ne.mype) then
              neighr = neighr + 1
              if (gr_msgbuffer) then
                 rprocs(neighr) = neigh(2,j,isg)
                 rtags (neighr) = neigh(1,j,isg)
              else
                 call MPI_IRECV(neigh_type(j,isg), & 
                      &                 1, & 
                      &                 MPI_INTEGER, & 
                      &                 neigh(2,j,isg), & 
                      &                 neigh(1,j,isg), & 
                      &                 MPI_COMM_WORLD, & 
                      &                 reqr(neighr), & 
                      &                 ierr)
              endif
           else
              neigh_type(j,isg) = nodetype(neigh(1,j,isg))
           end if
        else
           neigh_type(j,isg) = -1
        end if
     end do
  end do
  
! send messages if neighbor is off processor

  neighs = 0
  do isg = 1,lnblocks
     do jsend = 1,nfaces
        if(neigh(1,jsend,isg).gt.-1) then
           if(neigh(2,jsend,isg).ne.mype) then
              neighs = neighs + 1
              if (gr_msgbuffer) then 
                 sprocs(neighs) = neigh(2,jsend,isg)
                 stags (neighs) = isg
                 svals (neighs) = nodetype(isg)
              else
                 call MPI_SSEND(nodetype(isg), & 
                      &                 1, & 
                      &                 MPI_INTEGER, & 
                      &                 neigh(2,jsend,isg), &  ! PE TO SEND TO
                      &                 isg,         &  ! THIS IS THE TAG
                      &                 MPI_COMM_WORLD, & 
                      &                 ierr)
              endif
           end if
        end if
     end do
  end do
  
  if (.not.(gr_msgbuffer)) then
     if (neighr.gt.0) then
        call MPI_WAITALL (neighr, reqr, statr, ierr)
     end if
  else  
     
     call b_int_sendrcv(1025,1,neighs,sprocs,stags,svals, & 
          &                              neighr,rprocs,rtags,rvals)
     i = 0
     do isg = 1,lnblocks
        do j = 1,nfaces
           if(neigh(1,j,isg).gt.-1) then
              if(neigh(2,j,isg).ne.mype) then
                 i = i + 1
#ifdef ERRROR_CHECK
                 if ((rprocs(i).ne.neigh(2,j,isg)).or. & 
                      &                (rtags(i).ne.neigh(1,j,isg))) then
                    print *,mype,'Something horrible has happened.'
                 endif
#endif  
                 neigh_type(j,isg) = rvals(i)
              end if
           end if
        end do
     end do
     
  endif
  
  ! GET CHILD NODETYPES
  
  neighr = 0
  do isg = 1,lnblocks
     do j = 1,nchild
        if(child(1,j,isg).gt.-1) then
           if(child(2,j,isg).ne.mype) then
              neighr = neighr + 1
              if (gr_msgbuffer) then
                 rprocs(neighr) = child(2,j,isg)
                 rtags (neighr) = child(1,j,isg)
              else 
                 call MPI_IRECV(child_type(j,isg), & 
                      &                 1, & 
                      &                 MPI_INTEGER, & 
                      &                 child(2,j,isg), & 
                      &                 child(1,j,isg), & 
                      &                 MPI_COMM_WORLD, & 
                      &                 reqr(neighr), & 
                      &                 ierr)
              endif
           else
              child_type(j,isg) = nodetype(child(1,j,isg))
           end if
        else
           child_type(j,isg) = -1
        end if
     end do
  end do
  
  ! send messages if child is off processor
  
  neighs = 0
  do isg = 1,lnblocks
     if(parent(1,isg).gt.-1) then
        if(parent(2,isg).ne.mype) then
           neighs = neighs + 1
           if (gr_msgbuffer) then
              sprocs(neighs) = parent(2,isg)
              stags (neighs) = isg
              svals (neighs) = nodetype(isg)
           else
              call MPI_SSEND(nodetype(isg), & 
                   &               1, & 
                   &               MPI_INTEGER, & 
                   &               parent(2,isg), &  ! PE TO SEND TO
                   &               isg,           &  ! THIS IS THE TAG
                   &               MPI_COMM_WORLD, & 
                   &               ierr)
           endif
        end if
     end if
  end do
  
  if (.not.(gr_msgbuffer)) then
     if (neighr.gt.0) then
        call MPI_WAITALL (neighr, reqr, statr, ierr)
     end if
  else
     call b_int_sendrcv(1025,1,neighs,sprocs,stags,svals, & 
          &                              neighr,rprocs,rtags,rvals)
     i = 0
     do isg = 1,lnblocks
        do j = 1,nchild
           if(child(1,j,isg).gt.-1) then
              if(child(2,j,isg).ne.mype) then
                 i = i + 1
#ifdef ERRROR_CHECK
                 if ((rprocs(i).ne.child(2,j,isg)).or. & 
                      &                (rtags(i).ne.child(1,j,isg))) then
                    print *,mype,'Something horrible has happened.'
                 endif
#endif
                 child_type(j,isg) = rvals(i)
              end if
           end if
        end do
     end do
  end if
  
  return
end subroutine get_tree_nodetypes



