!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhExchangeTrees
!!
!! NAME
!!
!!  gr_bhExchangeTrees
!!
!!
!! SYNOPSIS
!!
!!   call gr_bhExchangeTrees()
!!
!! DESCRIPTION
!!   Determines levels up to which individual block trees need to be sent to
!!   different CPUs. For a given CPU, copies all trees up to appropriate levels
!!   to a single message and sends it to the CPU.
!!   
!!
!! ARGUMENTS
!!
!!
!!***



subroutine gr_bhExchangeTrees()

  use Grid_interface, ONLY : Grid_getLocalNumBlks, &
      Grid_getListOfBlocks, Grid_updateRefinement, &
      Grid_getBlkBoundBox
  use gr_bhData, ONLY : GR_TREE_BND_PERIODIC, &
    gr_bhLx, gr_bhLy, gr_bhLz, gr_bhLxHalf, gr_bhLyHalf, gr_bhLzHalf, &
    p_message, gr_bhTreeLevels, &
    gr_bhLocSentTreeLevels, gr_bhTreeNodetype, gr_bhTreeBCen, &
    gr_bhTreeArray, gr_bhComm, gr_bhTreeBS, gr_bhTreeLimAngle2, &
    gr_bhTreeMyPE, gr_bhTreeNumProcs, gr_bhBndType, &
    gr_bhTreeLrefine, gr_bhTreeCellSize, gr_bhTreeLnblocks, &
    gr_bhTreeNodeSize2, gr_bhLocRecvTreeLevels, gr_bhBlockTreePos, &
    GR_TREE_NSIZE, GR_TREE_IX, GR_TREE_IY, GR_TREE_IZ, &
    GR_TREE_BTP_POS ,GR_TREE_BTP_LEV, GR_TREE_BTP_C1, GR_TREE_BTP_C8, &
    gr_bhPhysMACComm_step, gr_bhBlockTreeNodeCen
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_bhLocalInterface, ONLY : gr_bhGetTreePos, gr_bhGetTreeSize, gr_bhMAC, &
    gr_bhPeriodicDr
  use tree, ONLY : nodetype, lrefine, child, mchild, lnblocks
  !use Logfile_interface, ONLY : Logfile_stamp

      
  implicit none
#include "constants.h"
#include "Flash_mpi.h"

  integer :: blockID, blockCount
  integer :: blockList(MAXBLOCKS)
  integer :: ierr, istat, i, j, k, ii, sgn
  integer :: lBID, rBID, l, ts, sp, pos, node_btp
  integer :: ms_size(0:gr_bhTreeNumProcs), mr_size(0:gr_bhTreeNumProcs)
  integer :: stats(MPI_STATUS_SIZE,2*(gr_bhTreeNumProcs-1)), reqs(2*(gr_bhTreeNumProcs-1))
  integer :: lcor, rcor, ls1, ls2, ls3, rs1, rs2, rs3
  integer :: multi(1:gr_bhTreeLevels), level
  integer, dimension(2,MDIM) :: nullBlkLimits = 0
  integer :: stack(1:(8**(gr_bhTreeLevels+1)-1)/7)
  integer, dimension(MDIM) :: nullpoint = 0
  type(p_message), save, allocatable :: messages_send(:), messages_recv(:)
  real ::  dist2
  real, allocatable  :: node(:)
  real, dimension(MDIM+2) :: dr ! two additional elements for dr^2 and 1/|dr|
  real, dimension(MDIM) :: drDiff, drGC
  logical :: node_accepted
  real, DIMENSION(:,:,:,:), POINTER :: solnData


  call Timers_start("exchange_tree")
  allocate(node(GR_TREE_NSIZE),  stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate node in gr_bhExchangeTrees.F90")

  ! determine levels of trees to be sent to different CPUs
  ! find distance of block centers minus both block space diagonals
  call Timers_start("et: det_gr_bhTreeLevels")
  do i = 0,gr_bhTreeNumProcs-1
    ms_size(i) = 0
    ! go through all local blocks
    do lBID = 1, gr_bhTreeLnblocks(gr_bhTreeMyPE)
      ! in this moment ignore parent blocks
      ! in future versions, if parent block would fullfill MAC, its local blocks will be not sent
      if (nodetype(lBID) .ne. 1) cycle

      gr_bhLocSentTreeLevels(lBID, i) = 0
      if (i /= gr_bhTreeMyPE) then

        ! TREE WALK: traverse the block_tree for the local block and evaluate
        ! MAC for each node and each remote block on a given CPU
        ! put root node on stack
        stack(1) = 1
        sp = 1
        do 
          ! take the node on the bottom of the stack
          node_btp = stack(sp)
          sp = sp - 1

          ! get the position in gr_bhTreeArray and the node level
          pos = gr_bhBlockTreePos(GR_TREE_BTP_POS, node_btp)
          level = gr_bhBlockTreePos(GR_TREE_BTP_LEV, node_btp)

          ! searching for the maximum level necessary 
          if (level > gr_bhLocSentTreeLevels(lBID, i)) then
            gr_bhLocSentTreeLevels(lBID, i) = level
          endif

          ! reached level of individual cells: exit the tree walk
          if (level == gr_bhTreeLevels) goto 666

          ! construct the node array
          node = gr_bhTreeArray(gr_bhTreeMyPE, lBID)%p(pos:pos+GR_TREE_NSIZE-1)

          ! find remote block with minimum distance from the node
          do rBID = 1, gr_bhTreeLnblocks(i) ! remote blocks (on CPU i)
            if (gr_bhTreeNodetype(rBID, i) .ne. 1) cycle

            ! compute Bmin vector: distance between node mass centre and 
            ! the closest point from the remote block
            ! 1. local_nodeMass_centre - remote_block_centre)
            dr(IAXIS) = gr_bhTreeArray(gr_bhTreeMyPE, lBID)%p(pos+GR_TREE_IX-1) &
            &         - gr_bhTreeBCen(IAXIS,rBID,i)
            dr(JAXIS) = gr_bhTreeArray(gr_bhTreeMyPE, lBID)%p(pos+GR_TREE_IY-1) &
            &         - gr_bhTreeBCen(JAXIS,rBID,i)
            dr(KAXIS) = gr_bhTreeArray(gr_bhTreeMyPE, lBID)%p(pos+GR_TREE_IZ-1) &
            &         - gr_bhTreeBCen(KAXIS,rBID,i)
            ! 2. find the minimum distance among copies of Ewald cells
            call gr_bhPeriodicDr(dr)
            ! 3. subtract sgn(dr)*0.5*remote_block_size, make sure dr does not change sign
            if (dr(IAXIS) < 0) then
              dr(IAXIS) = min(0.0, dr(IAXIS) + 0.5 * gr_bhTreeBS &
              & * gr_bhTreeCellSize(gr_bhTreeLrefine(rBID, i), IAXIS))
            else
              dr(IAXIS) = max(0.0, dr(IAXIS) - 0.5 * gr_bhTreeBS &
              & * gr_bhTreeCellSize(gr_bhTreeLrefine(rBID, i), IAXIS))
            endif
            if (dr(JAXIS) < 0) then
              dr(JAXIS) = min(0.0, dr(JAXIS) + 0.5 * gr_bhTreeBS &
              & * gr_bhTreeCellSize(gr_bhTreeLrefine(rBID, i), JAXIS))
            else
              dr(JAXIS) = max(0.0, dr(JAXIS) - 0.5 * gr_bhTreeBS &
              & * gr_bhTreeCellSize(gr_bhTreeLrefine(rBID, i), JAXIS))
            endif
            if (dr(KAXIS) < 0) then
              dr(KAXIS) = min(0.0, dr(KAXIS) + 0.5 * gr_bhTreeBS &
              & * gr_bhTreeCellSize(gr_bhTreeLrefine(rBID, i), KAXIS))
            else
              dr(KAXIS) = max(0.0, dr(KAXIS) - 0.5 * gr_bhTreeBS &
              & * gr_bhTreeCellSize(gr_bhTreeLrefine(rBID, i), KAXIS))
            endif

            dr(MDIM+1)  = dr(IAXIS)*dr(IAXIS) + dr(JAXIS)*dr(JAXIS) + dr(KAXIS)*dr(KAXIS)
            dr(MDIM+2) = sqrt(1.0 / (dr(MDIM+1) + 1D-99))

            ! position vector to the node geometric centre
            ! drDiff is a vector from the node mass centre to its geometrical centre
            drDiff(:) = gr_bhTreeBCen(:, lBID, gr_bhTreeMyPE) &
            & + gr_bhTreeCellSize(gr_bhTreeLrefine(lBID, gr_bhTreeMyPE), :) &
            & * gr_bhBlockTreeNodeCen(:, node_btp) &
            & - node(GR_TREE_IX:GR_TREE_IZ)
            drGC = dr(1:MDIM) + drDiff

            ! test the multipole acceptence criterion
            node_accepted = gr_bhMAC(gr_bhPhysMACComm_step, node &
            & , gr_bhTreeNodeSize2(level+gr_bhTreeLrefine(lBID,gr_bhTreeMyPE)) &
            & , drGC, dr, -1, nullpoint, nullBlkLimits, solnData)

            if (.not. node_accepted) exit

          enddo


   


          if (.not. node_accepted) then
            ! put all children on the stack
            do l = GR_TREE_BTP_C1,GR_TREE_BTP_C8
              sp = sp + 1
              stack(sp) = gr_bhBlockTreePos(l, node_btp)
            enddo
          endif
          
          ! stack is empty - exiting
          if (sp == 0) exit

        enddo ! END OF TREE WALK

666     continue
        ! add tree size to the ms_size
        ms_size(i) = ms_size(i) + gr_bhGetTreeSize(gr_bhLocSentTreeLevels(lBID, i))
      endif
    enddo
  enddo
  call Timers_stop("et: det_gr_bhTreeLevels")
  
  call Timers_start("et: comm_gr_bhTreeLevels")

  !if (gr_bhTreeMyPE == MASTER_PE) then
  !   call Logfile_stamp( "ET before STL", "[BHTree]")
  !end if
  ! sent sentTreeLevels to recvTreeLevels :)
  j = 1 ! tag of the message, use also for indexing reqs and stats arrays
  do i = 0,gr_bhTreeNumProcs-1
    if (i /= gr_bhTreeMyPE) then
      call MPI_Irecv(gr_bhLocRecvTreeLevels(1,i), MAXBLOCKS, FLASH_INTEGER, i, 1, gr_bhComm, reqs(j), ierr)
      call MPI_Isend(gr_bhLocSentTreeLevels(1,i), MAXBLOCKS, FLASH_INTEGER, i, 1, gr_bhComm, reqs(gr_bhTreeNumProcs-1+j), ierr)
      j = j + 1
    endif
  enddo
  call MPI_WaitAll(2*(gr_bhTreeNumProcs-1), reqs, stats, ierr)
  

  !if (gr_bhTreeMyPE == MASTER_PE) then
  !   call Logfile_stamp( "ET after STL", "[BHTree]")
  !end if

  call Timers_stop("et: comm_gr_bhTreeLevels")
  call Timers_start("et: comm_trees")

  ! build messages to be sent, fill them with local trees
  allocate(messages_send(0:gr_bhTreeNumProcs), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate messages_send in gr_bhExchangeTrees.F90")
  do i = 0,gr_bhTreeNumProcs-1
    if (i /= gr_bhTreeMyPE) then
      ! create message to be sent to i-th CPU
      nullify(messages_send(i)%p)
      allocate(messages_send(i)%p(1:ms_size(i)), stat=istat)
      if (istat /= 0) call Driver_abortFlash("could not allocate messages_send in gr_bhExchangeTrees.F90")
      j = 1
      do lBID = 1, gr_bhTreeLnblocks(gr_bhTreeMyPE)
        if (nodetype(lBID) .ne. 1) cycle
        do k = 1, gr_bhGetTreeSize(gr_bhLocSentTreeLevels(lBID, i))
          messages_send(i)%p(j) = gr_bhTreeArray(gr_bhTreeMyPE, lBID)%p(k)
          j = j + 1
        enddo
      enddo
    endif
  enddo

  
  ! prepare space for messages to be recieved
  allocate(messages_recv(0:gr_bhTreeNumProcs), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate messages_recv in gr_bhExchangeTrees.F90")
  do i = 0,gr_bhTreeNumProcs-1
    if (i /= gr_bhTreeMyPE) then
      j = 0
      do rBID = 1, gr_bhTreeLnblocks(i)
        if (gr_bhTreeNodetype(rBID, i) .ne. 1) cycle
        j = j + gr_bhGetTreeSize(gr_bhLocRecvTreeLevels(rBID, i))
      enddo
      mr_size(i) = j
      nullify(messages_recv(i)%p)
      allocate(messages_recv(i)%p(1:mr_size(i)))
      if (istat /= 0) call Driver_abortFlash("could not allocate messages_recv in gr_bhExchangeTrees.F90")
    endif
  enddo

  ! to ensure buffers for messages are allocated
  ! may not be necessary
  call MPI_Barrier(gr_bhComm, ierr)
  !if (gr_bhTreeMyPE == MASTER_PE) then
  !   call Logfile_stamp( "ET tree mes allocated", "[BHTree]")
  !end if
  
  ! exchange messages
  j = 1 ! tag of the message, use also for indexing reqs and stats arrays
  do i = 0,gr_bhTreeNumProcs-1
    if (i /= gr_bhTreeMyPE) then
      call MPI_Irecv(messages_recv(i)%p, mr_size(i), FLASH_REAL, i, 1, gr_bhComm, reqs(j), ierr)
      call MPI_Isend(messages_send(i)%p, ms_size(i), FLASH_REAL, i, 1, gr_bhComm, reqs(gr_bhTreeNumProcs-1+j), ierr)
      j = j + 1
    endif
  enddo
  call MPI_WaitAll(2*(gr_bhTreeNumProcs-1), reqs, stats, ierr)
  
  !if (gr_bhTreeMyPE == MASTER_PE) then
  !   call Logfile_stamp( "ET tree exchanged", "[BHTree]")
  !end if
  
  ! copy data from messages to trees
  do i = 0,gr_bhTreeNumProcs-1
    if (i /= gr_bhTreeMyPE) then
      j = 1
      do rBID = 1, gr_bhTreeLnblocks(i)
        if (gr_bhTreeNodetype(rBID, i) .ne. 1) cycle
        ts = gr_bhGetTreeSize(gr_bhLocRecvTreeLevels(rBID, i))
        nullify(gr_bhTreeArray(i,rBID)%p)
        allocate(gr_bhTreeArray(i,rBID)%p(1:ts), stat=istat)
        if (istat /= 0) call Driver_abortFlash("could not allocate tree in gr_bhExchangeTrees.F90")
        do k = 1, ts
          gr_bhTreeArray(i, rBID)%p(k) = messages_recv(i)%p(j)
          j = j + 1
        enddo
      enddo
    endif
  enddo

  ! destroy messages
  do i = 0,gr_bhTreeNumProcs-1
    if (i /= gr_bhTreeMyPE) then
      deallocate(messages_send(i)%p)
      deallocate(messages_recv(i)%p)
    endif
  enddo
  deallocate(messages_send)
  deallocate(messages_recv)
  call Timers_stop("et: comm_trees")
  
  call Timers_stop("exchange_tree")
  deallocate(node)

  return
end subroutine gr_bhExchangeTrees


