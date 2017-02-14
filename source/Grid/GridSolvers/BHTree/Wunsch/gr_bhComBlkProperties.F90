!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhComBlkProperties
!!
!! NAME
!!
!!  gr_bhComBlkProperties
!!
!!
!! SYNOPSIS
!!
!!   call gr_bhComBlkProperties()
!!
!! DESCRIPTION
!!
!!   Communicates global block properties: nodetype, lrefine, child,
!!     and position of the block centre (gr_bhTreeBCen).
!!
!! ARGUMENTS
!!
!!             
!!
!!***



subroutine gr_bhComBlkProperties()

  use tree, ONLY : nodetype, maxblocks_tr, mchild, mfaces &
                &, lrefine, child, neigh, lnblocks
  use Grid_interface, ONLY : Grid_getCellCoords, &
      Grid_getListOfBlocks
  use gr_bhData, ONLY : gr_bhTreeNodetype, gr_bhTreeLrefine, &
    gr_bhTreeNumProcs, gr_bhTreeMaxlnblocks, gr_bhComm, gr_bhTreeLrefine, &
    gr_bhLocCoords, gr_bhTreeChild, gr_bhTreeNeigh, gr_bhTreeBCen, &
    gr_bhTreeFirstLevBlocks, gr_bhTreeLnblocks, gr_bhTreeBlockList, &
    gr_bhTreeNFirstLev, gr_bhTreeBS, gr_bhTreeMyPE
  use gr_bhLocalInterface, ONLY : gr_bhGetTreeSize
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash
  use Timers_interface, ONLY : Timers_start, Timers_stop
  !use Logfile_interface, ONLY : Logfile_stamp
      
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer :: ierr, istat, i, j, k, l, jjj
  integer :: block, int_mes_size
  integer :: stats(MPI_STATUS_SIZE,4*(gr_bhTreeNumProcs-1)), reqs(4*(gr_bhTreeNumProcs-1))
  integer,allocatable :: sent_int_mes(:), recv_int_mes(:,:)
  real,allocatable :: sent_real_mes(:), recv_real_mes(:,:)
  integer ii, blockCount
  logical       :: gcell = .false.
  integer :: blockList(MAXBLOCKS)
  character(len=256) :: strBuff

  call Timers_start("comm_blk_props")

  ! set maximum amount of blocks on 1 CPU
  call MPI_ALLREDUCE (lnblocks,gr_bhTreeMaxlnblocks,1,FLASH_INTEGER, & 
     &                    MPI_MAX,gr_bhComm,ierr)

  do i = 0, gr_bhTreeNumProcs-1
    do j = 1, maxblocks_tr
        gr_bhTreeNodetype(j,i) = -1
        gr_bhTreeLrefine(j,i) = -1
      enddo
  enddo
  
  
  ! get all local blocks coordinates
  call Grid_getListOfBlocks(ALL_BLKS,blockList,blockCount) ! get list of LEAF blocks
  do i = 1, blockCount
    block = blockList(i)

    ! get information about coordinates of the block
    call Grid_getCellCoords(IAXIS, block, CENTER, gcell, gr_bhLocCoords(:, IAXIS, block), gr_bhTreeBS)
    call Grid_getCellCoords(JAXIS, block, CENTER, gcell, gr_bhLocCoords(:, JAXIS, block), gr_bhTreeBS)
    call Grid_getCellCoords(KAXIS, block, CENTER, gcell, gr_bhLocCoords(:, KAXIS, block), gr_bhTreeBS)

    gr_bhLocCoords(gr_bhTreeBS+1, IAXIS, i) = 0.5*(gr_bhLocCoords(1, IAXIS, block)+gr_bhLocCoords(gr_bhTreeBS, IAXIS, block))
    gr_bhLocCoords(gr_bhTreeBS+1, JAXIS, i) = 0.5*(gr_bhLocCoords(1, JAXIS, block)+gr_bhLocCoords(gr_bhTreeBS, JAXIS, block))
    gr_bhLocCoords(gr_bhTreeBS+1, KAXIS, i) = 0.5*(gr_bhLocCoords(1, KAXIS, block)+gr_bhLocCoords(gr_bhTreeBS, KAXIS, block))
  enddo

  ! allocate messages
  int_mes_size = (3+2*MCHILD)*gr_bhTreeMaxlnblocks+1
  allocate(sent_int_mes(1:int_mes_size), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate sent_int_mes in treeComBlkProperties")
  allocate(recv_int_mes(1:int_mes_size, 0:gr_bhTreeNumProcs-1), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate recv_int_mes in treeComBlkProperties")
  allocate(sent_real_mes(1:MDIM*gr_bhTreeMaxlnblocks), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate sent_real_mes in treeComBlkProperties")
  allocate(recv_real_mes(1:MDIM*gr_bhTreeMaxlnblocks, 0:gr_bhTreeNumProcs-1), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate recv_real_mes in treeComBlkProperties")

  ! copy data messages:
  ! integer message: indeces(gr_bhTreeMaxlnblocks), nodetype(gr_bhTreeMaxlnblocks)
  !                , lrefine(gr_bhTreeMaxlnblocks), child(2, MCHILD, gr_bhTreeMaxlnblocks)
  ! real message: gr_bhLocCoords(NXB+1,MDIM,gr_bhTreeMaxlnblocks)
  j = 1
  do i = 1, lnblocks
    if (nodetype(i) .ne. -1) then
      sent_int_mes(j) = i
      j = j + 1
    endif
  enddo
  do i = j,gr_bhTreeMaxlnblocks ! fill the rest of indeces with -1
    sent_int_mes(i) = -1
  enddo
  do i = 1,gr_bhTreeMaxlnblocks ! copy nodetype and lrefine into sent_int_mes
    if (sent_int_mes(i) .ne. -1) then
      sent_int_mes(gr_bhTreeMaxlnblocks+i) = nodetype(sent_int_mes(i))
      sent_int_mes(2*gr_bhTreeMaxlnblocks+i) = lrefine(sent_int_mes(i))
      do j = 1, mchild
        sent_int_mes(3*gr_bhTreeMaxlnblocks + (i-1)*2*mchild + 2*(j-1) + 1) = child(1,j,sent_int_mes(i))
        sent_int_mes(3*gr_bhTreeMaxlnblocks + (i-1)*2*mchild + 2*(j-1) + 2) = child(2,j,sent_int_mes(i))
      enddo
    endif
  enddo
  sent_int_mes(int_mes_size) = lnblocks

  ! real message: (gr_bhTreeMaxlnblocks), nodetype(gr_bhTreeMaxlnblocks), lrefine(gr_bhTreeMaxlnblocks)
  do i = 1,gr_bhTreeMaxlnblocks ! copy nodetype and lrefine into sent_int_mes
    if (sent_int_mes(i) .ne. -1) then
      do j = 1,MDIM
          sent_real_mes((i-1)*MDIM+j) = gr_bhLocCoords(NXB+1,j,sent_int_mes(i))
      enddo
    endif
  enddo


  call MPI_Barrier(gr_bhComm, ierr) ! to ensure that messages are allocated before received

  ! send messages
  j = 1 ! tag of the message, use also for indexing reqs and stats arrays
  do i = 0,gr_bhTreeNumProcs-1
    if (i /= gr_bhTreeMyPE) then
      call MPI_Irecv(recv_int_mes(1,i), int_mes_size, FLASH_INTEGER, i, 1, gr_bhComm, reqs(j), ierr)
      call MPI_Isend(sent_int_mes(1), int_mes_size, FLASH_INTEGER, i, 1, gr_bhComm, reqs(2*(gr_bhTreeNumProcs-1)+j), ierr)
      j = j + 1
      call MPI_Irecv(recv_real_mes(1,i), MDIM*gr_bhTreeMaxlnblocks, FLASH_REAL, i, 1, gr_bhComm, reqs(j), ierr)
      call MPI_Isend(sent_real_mes(1), MDIM*gr_bhTreeMaxlnblocks, FLASH_REAL, i, 1, gr_bhComm, &
           reqs(2*(gr_bhTreeNumProcs-1)+j), ierr)
      j = j + 1
    else
      recv_int_mes(:,i) = sent_int_mes(:)
      recv_real_mes(:,i) = sent_real_mes(:)
    endif
  enddo
  call MPI_WaitAll(4*(gr_bhTreeNumProcs-1), reqs, stats, ierr)

  ! copy data back from messages to global gr_bh* structures
  do l = 0,gr_bhTreeNumProcs-1
    do i = 1,gr_bhTreeMaxlnblocks
      if (recv_int_mes(i,l) .ne. -1) then
        gr_bhTreeNodetype(recv_int_mes(i,l),l) = recv_int_mes(gr_bhTreeMaxlnblocks+i,l)
        gr_bhTreeLrefine(recv_int_mes(i,l),l)  = recv_int_mes(2*gr_bhTreeMaxlnblocks+i,l)
        
        do j = 1, mchild
          gr_bhTreeChild(1,j,recv_int_mes(i,l),l) = recv_int_mes(3*gr_bhTreeMaxlnblocks + (i-1)*2*mchild + 2*(j-1) + 1,l)
          gr_bhTreeChild(2,j,recv_int_mes(i,l),l) = recv_int_mes(3*gr_bhTreeMaxlnblocks + (i-1)*2*mchild + 2*(j-1) + 2,l)
        enddo
      endif
    enddo
    gr_bhTreeLnblocks(l) = recv_int_mes(int_mes_size,l)
  enddo

  do l = 0,gr_bhTreeNumProcs-1
    do i = 1,gr_bhTreeMaxlnblocks
      if (recv_int_mes(i,l) .ne. -1) then
        do j = 1,MDIM
            gr_bhTreeBCen(j,recv_int_mes(i,l),l) = recv_real_mes((i-1)*MDIM+j,l)
        enddo
      endif
    enddo
  enddo



  call MPI_Barrier(gr_bhComm, ierr) ! to ensure that messages are allocated before received

  deallocate(gr_bhTreeBlocklist)
  allocate(gr_bhTreeBlocklist(1:gr_bhTreeMaxlnblocks, 0:gr_bhTreeNumProcs-1), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhTreeBlocklist in treeComBlkProperties")
  do l = 0,gr_bhTreeNumProcs-1
    gr_bhTreeBlocklist(:,l) = recv_int_mes(1:gr_bhTreeMaxlnblocks,l)
  enddo
  
  ! and destroy messages
  deallocate(sent_int_mes)
  deallocate(recv_int_mes)
  deallocate(sent_real_mes)
  deallocate(recv_real_mes)


  ! find blocks on the top level
  k = 0
  do i = 0,gr_bhTreeNumProcs-1
    do j = 1, gr_bhTreeLnblocks(i)
      if (gr_bhTreeLrefine(j, i) .eq. 1) then
        k = k + 1
        if (k > gr_bhTreeNFirstLev) then
          call Driver_abortFlash("gr_bhTreeNFirstLev exceeded in gr_bhComBlkProperties.F90")
        endif
        gr_bhTreeFirstLevBlocks(1, k) = j
        gr_bhTreeFirstLevBlocks(2, k) = i
      endif
    enddo
  enddo


  call Timers_stop("comm_blk_props")

  !if (gr_bhTreeMyPE == MASTER_PE) then
  !   call Logfile_stamp("block properties exchanged", "[BHTree]")
  !end if

  return
end subroutine gr_bhComBlkProperties



