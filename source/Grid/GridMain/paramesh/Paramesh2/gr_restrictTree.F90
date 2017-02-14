!!****if* source/Grid/GridMain/paramesh/Paramesh2/gr_restrictTree
!!
!! NAME
!!  gr_restrictTree
!!
!! SYNOPSIS
!!
!!  gr_restrictTree()
!!
!!  
!! DESCRIPTION 
!!
!!  This subroutine restrict the tree all the way to the coarsest
!!  refinement level. It is usually invoked before writing data to
!!  and output file. This is mostly for visualization purposes to
!!  be able to look at different levels of resolution
!!
!!
!! ARGUMENTS   
!!
!!***


subroutine gr_restrictTree()
  
  use physicaldata
  use tree
  use workspace
  implicit none
  include 'Flash_mpi.h'
  
  
  integer nodetype_old(maxblocks_tr),level_at,lmax,lmax2
  integer nprocs,mype
  integer i,ierr
  logical :: first_call=.true.

  if(first_call) then
     first_call = .false.
     return
  end if
  
  lmax = -1
  do i = 1,lnblocks
     if (nodetype(i).eq.1) then
        lmax = max(lrefine(i),lmax)
     end if
  end do
  
  call MPI_COMM_SIZE (MPI_COMM_WORLD,nprocs,ierr)
  call MPI_COMM_RANK (MPI_COMM_WORLD,mype,  ierr)
  call MPI_ALLREDUCE (lmax,lmax2,1,MPI_INTEGER, & 
       &     MPI_MAX,MPI_COMM_WORLD,ierr)
  lmax = lmax2
  
  nodetype_old(:) = nodetype(:)
  
  do level_at = lmax,2,-1
     
     call amr_restrict(mype,1,0)
     
     do i = 1, lnblocks
        if (nodetype(i).eq.1.and.lrefine(i).eq.level_at) & 
             &           nodetype(i) = -1
        if (nodetype(i).eq.2.and.lrefine(i).eq.level_at-1) & 
             &           nodetype(i) = 1
        if (nodetype(i).eq.3.and.lrefine(i).eq.level_at-2) & 
             &           nodetype(i) = 2
     end do
     
     do i = 1,lnblocks
        if (nodetype(i).eq.100) nodetype(i) = 1
     end do
     
     ! reset neighbor nodetypes and child nodetypes
     call get_tree_nodetypes (nprocs,mype)
     
  end do
  
  nodetype(:) = nodetype_old(:)
  ! reset neighbor nodetypes and child nodetypes
  call get_tree_nodetypes (nprocs,mype)
  
  return
end subroutine gr_restrictTree

