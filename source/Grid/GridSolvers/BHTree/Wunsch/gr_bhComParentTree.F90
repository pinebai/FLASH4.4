!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhComParentTree
!!
!! NAME
!!
!!  gr_bhComParentTree
!!
!!
!! SYNOPSIS
!!
!!   gr_bhComParentTree(
!!        integer, intent(IN) :: level
!!        )
!!
!! DESCRIPTION
!!
!!   Communicates the array gr_bhTreeParentTree which includes masses and mass 
!!     centres of all blocks. Values for blocks at a given refinement level are
!!     distributed to all CPUs. If level = -1, values for all LEAF blocks
!!     are distributed.
!!
!! ARGUMENTS
!!
!!  level : this refinement level of the parent tree is communicated
!!          if level == -1: all leaf blocks are communicated
!!
!!***



subroutine gr_bhComParentTree(level)

  use gr_bhData, ONLY : gr_bhTreeNumProcs, gr_bhTreeLnblocks, &
    gr_bhTreeNodetype, gr_bhTreeLrefine, gr_bhComm, gr_bhTreeLnblocks, &
    gr_bhTreeParentTree, gr_bhTreeMyPE, GR_TREE_NSIZE
  use gr_bhLocalInterface, ONLY : gr_bhGetTreeSize
      
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer, intent(IN) :: level
  integer :: istat, i, j, l
  integer :: cpu_size, mes_size
  integer :: stats(MPI_STATUS_SIZE,2*(gr_bhTreeNumProcs-1)), reqs(2*(gr_bhTreeNumProcs-1))
  real, allocatable :: sent_mes(:,:), recv_mes(:,:,:)
  integer ii

  ! determine size of messages
  mes_size = 0
  do l = 0, gr_bhTreeNumProcs-1
    cpu_size = 0
    do i = 1, gr_bhTreeLnblocks(l)
      if (level .eq. -1) then
        ! exchanging leaf blocks
        if (gr_bhTreeNodetype(i,l) .eq. 1) cpu_size = cpu_size + 1
      else
        ! exchanging all blocks at a given level
        if (gr_bhTreeLrefine(i,l) .eq. level) cpu_size = cpu_size + 1
      endif
    enddo
    if (cpu_size .gt. mes_size) mes_size = cpu_size
  enddo
  if (mes_size == 0) mes_size = 1

  ! allocate arrays for messages
  allocate(sent_mes(1:GR_TREE_NSIZE,1:mes_size), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate sent_real_mes in treeComParentTree")
  allocate(recv_mes(1:GR_TREE_NSIZE,1:mes_size, 0:gr_bhTreeNumProcs-1), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate recv_real_mes in treeComParentTree")

  ! coppy data from gr_bhLocParentTree into messages
  call MPI_Barrier(gr_bhComm, istat) ! to ensure all recv_* buffers are allocated
  j = 1
  do i = 1,gr_bhTreeLnblocks(gr_bhTreeMyPE)
    if (level .eq. -1) then
      ! exchanging leaf blocks
      if (gr_bhTreeNodetype(i,gr_bhTreeMyPE) .eq. 1) then
        sent_mes(:,j) = gr_bhTreeParentTree(:, i, gr_bhTreeMyPE)
        j = j + 1
      endif
    else
      ! exchanging all blocks at a given level
      if (gr_bhTreeLrefine(i,gr_bhTreeMyPE) .eq. level) then
        sent_mes(:,j) = gr_bhTreeParentTree(:, i, gr_bhTreeMyPE)
        j = j + 1
      endif
    endif
  enddo
  
  ! send messages
  j = 1 ! tag of the message, use also for indexing reqs and stats arrays
  do i = 0,gr_bhTreeNumProcs-1
    if (i /= gr_bhTreeMyPE) then
      call MPI_Irecv(recv_mes(1,1,i), GR_TREE_NSIZE*mes_size, FLASH_REAL, i, 1, gr_bhComm, reqs(j), istat)
      call MPI_Isend(sent_mes(1,1), GR_TREE_NSIZE*mes_size, FLASH_REAL, i, 1, gr_bhComm, reqs((gr_bhTreeNumProcs-1)+j), istat)
      j = j + 1
    endif
  enddo
  call MPI_WaitAll(2*(gr_bhTreeNumProcs-1), reqs, stats, istat)
  
  ! copy data back from messages into gr_bhTreeParentTree
  do l = 0,gr_bhTreeNumProcs-1
    if (l /= gr_bhTreeMyPE) then
      j = 1
      do i = 1,gr_bhTreeLnblocks(l)
        if (level .eq. -1) then
          ! exchanging leaf blocks
          if (gr_bhTreeNodetype(i,l) .eq. 1) then
            gr_bhTreeParentTree(:, i, l) = recv_mes(:,j, l)
            j = j + 1
          endif
        else
          ! exchanging all blocks at a given level
          if (gr_bhTreeLrefine(i,l) .eq. level) then
            gr_bhTreeParentTree(:, i, l) = recv_mes(:,j, l)
            j = j + 1
          endif
        endif
      enddo
    endif
  enddo

  ! destroy messages
  deallocate(sent_mes)
  deallocate(recv_mes)
  call MPI_Barrier(gr_bhComm, istat) ! to ensure all Parent Trees are at place before the global tree is built

  return
end subroutine gr_bhComParentTree




