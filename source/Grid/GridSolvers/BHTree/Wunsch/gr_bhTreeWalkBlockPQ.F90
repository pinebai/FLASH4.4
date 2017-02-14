!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhTreeWalkBlockPQ
!!
!! NAME
!!
!!  gr_bhTreeWalkBlock
!!
!!
!! SYNOPSIS
!!
!!   gr_bhTreeWalkBlock(
!!          integer,intent(in) :: block
!!          )
!!
!! DESCRIPTION
!!
!!   Traverses the tree in a specific block. Tree is traversed using the
!!   priority-queue, nodes with the largest estimated error of the contribution
!!   are precessed first. It is used by the SumSquare MAC (see Salmon&Warren,
!!   1994, J.Comp.Phys, 136, 155 for details). This algorithm is not fully
!!   implemented in this version (specifically, communication in
!!   gr_bhExchangeTrees is not compatible with the included implementation)
!!
!! ARGUMENTS
!!
!!  block    - ID of a block where the contribution is calculated
!!
!!***



subroutine gr_bhTreeWalkBlockPQ(block)

!   use Grid_interface, ONLY : Grid_getLocalNumBlks, &
!     Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getDeltas
!   use gr_bhLocalInterface, ONLY : gr_bhLeafContrib, gr_bhParentContrib &
!     , gr_bhStartBlock, gr_bhFinalizeBlock, gr_bhMAC, gr_bhBotNodeContrib &
!     , gr_bhNodeContrib, gr_bhSelfContrib, gr_bhPQInsert, gr_bhPQPull
! 
!   use gr_bhData, ONLY : gr_bhLocCoords, gr_bhTreeMyPE, &
!     gr_bhTreeLrefine, gr_bhTreeNodeSize2, gr_bhPhysMACTW_step, &
!     gr_bhTreeDcount, gr_bhTreeNFirstLev, gr_bhTreeChild, gr_bhTreeNodetype, &
!     gr_bhTreeFirstLevBlocks, gr_bhTreeParentTree, &
!     gr_bhTreeCellSize, gr_bhTreeArray, &
!     GR_TREE_NSIZE, GR_TREE_BNSIZE, GR_TREE_IX, GR_TREE_IY, GR_TREE_IZ, &
!     GR_TREE_BTP_POS ,GR_TREE_BTP_LEV, GR_TREE_BTP_I, GR_TREE_BTP_J, GR_TREE_BTP_K, &
!     GR_TREE_BTP_C1, GR_TREE_BTP_C8, gr_bhLocRecvTreeLevels, &
!     gr_bhBlockTreePos, gr_bhTreeNodeSize, GR_TREE_INT_CC, GR_TREE_INT_CN, GR_TREE_INT_CP, &
!     GR_TREE_INT_SELF, gr_bhPQSize, gr_bhPQSumSquare, gr_bhPriorityQueue, gr_bhPQNull, & 
!     gr_bhTreeBCen, gr_bhTreeLevels, gr_bhTreeLoff, GR_TREE_IM
!   use tree, ONLY : nodetype, child, mchild
!   use Driver_interface, ONLY : Driver_abortFlash
  implicit none
#include "Flash.h"
#include "constants.h"
#include "Flash_mpi.h"
  integer, intent(IN) :: block
!   integer       :: i, j, k, l, cpu, tr, sp, istat, level
!   integer       :: ii, jj, kk, point_mbLp1(MDIM), ipq
!   integer       :: node_btp, pos, int_type
!   real          :: dcount(1:3), perr, ss
!   real, POINTER, DIMENSION(:,:,:,:) :: solnData
!   integer, dimension(2,MDIM)   :: blkLimits,blkLimitsGC
!   real, dimension(MDIM+2) :: dr ! two additional elements dr^2 and 1/|dr|
!   integer, dimension(MDIM) :: point
!   real, dimension(MDIM) :: phys_point, del
!   real, allocatable  :: node(:), botnode(:)
!   logical       :: node_accepted, gcell = .false.
!   integer :: lev
!   real :: rc, r2c, r3c, m, x0, y0, z0, b3, b2, dvol, b3_cell, b2_cell

  ! This feature is not implemented yet
  return

!  ! allocate the node array
!  allocate(node(GR_TREE_NSIZE), stat=istat)
!  if (istat /= 0) call Driver_abortFlash("could not allocate node in gr_bhTreeWalkBlockPQ.F90")
!  allocate(botnode(GR_TREE_BNSIZE), stat=istat)
!  if (istat /= 0) call Driver_abortFlash("could not allocate botnode in gr_bhTreeWalkBlockPQ.F90")
!
!  ! get information about the block
!  ! indeces necessary to write values into solnData
!  call Grid_getBlkIndexLimits(block,blkLimits,blkLimitsGC) 
!  ! pointer to the density and gp field
!  call Grid_getBlkPtr(block,solnData,CENTER) 
!
!  ! call physics routines to do block initialization
!  call gr_bhStartBlock(block, blkLimits, solnData)
!  
!  ! reset distance counter
!  dcount = 0.0
!  node = gr_bhTreeParentTree(:, block, gr_bhTreeMyPE)
!  !print *, "node = ", node
!  !print *, "BCen = ", gr_bhTreeBCen(:, block, gr_bhTreeMyPE)
!  x0 = node(GR_TREE_IX)
!  y0 = node(GR_TREE_IY)
!  z0 = node(GR_TREE_IZ)
!  b3 = 0.0
!  b2 = 0.0
!  call Grid_getDeltas(block, del)
!  do kk = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS) ! this block
!    point(KAXIS) = kk
!    point_mbLp1(KAXIS) = kk - blkLimits(LOW, KAXIS) + 1
!    phys_point(KAXIS) = gr_bhLocCoords(point_mbLp1(KAXIS), KAXIS, block)
!    do jj = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
!      point(JAXIS) = jj
!      point_mbLp1(JAXIS) = jj - blkLimits(LOW, JAXIS) + 1
!      phys_point(JAXIS) = gr_bhLocCoords(point_mbLp1(JAXIS), JAXIS, block)
!      do ii = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
!        point(IAXIS) = ii
!        point_mbLp1(IAXIS) = ii - blkLimits(LOW, IAXIS) + 1
!        phys_point(IAXIS) = gr_bhLocCoords(point_mbLp1(IAXIS), IAXIS, block)
!        call Grid_getSingleCellVol(block, INTERIOR, point, dvol)
!        m = solnData(DENS_VAR, ii, jj, kk)*dvol
!        r2c = (phys_point(IAXIS)-x0)**2 + (phys_point(JAXIS)-y0)**2 + (phys_point(KAXIS)-z0)**2
!        rc  = sqrt(r2c)
!        r3c = rc*r2c
!        b2_cell = (1.0/12.0)*m*(del(IAXIS)**2 + del(JAXIS)**2 + del(KAXIS)**2)
!        !b3_cell = GRAV_B3_CONST * m * del(IAXIS)**3
!        !b3 = b3 + m*abs(r3c) + b3_cell + 3*b2_cell*rc
!        b2 = b2 + m*r2c + b2_cell
!      enddo
!    enddo
!  enddo
!  !print *, "B3: ", b3, node(GRAV_IB3), sqrt(node(GRAV_IB2)**3/node(GR_TREE_IM))
!  !print *, "B2: ", b2, node(GRAV_IB2)
!
!  !do i = 0, gr_bhTreeLevels-1
!  !  do j = 0,8**i-1
!  !    pos = gr_bhTreeLoff(i)+GR_TREE_NSIZE*j
!  !    node = gr_bhTreeArray(gr_bhTreeMyPE, block)%p(pos:pos+GR_TREE_NSIZE-1)
!  !    print *, "lev,node,pos,m,B2,B3: ", i, j,pos,node(GR_TREE_IM), &
!  !    & node(GRAV_IB2), node(GRAV_IB3), sqrt(node(GRAV_IB2)**3/node(GR_TREE_IM))
!  !  enddo
!  !enddo
!
!  do kk = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS) ! this block
!    point(KAXIS) = kk
!    point_mbLp1(KAXIS) = kk - blkLimits(LOW, KAXIS) + 1
!    phys_point(KAXIS) = gr_bhLocCoords(point_mbLp1(KAXIS), KAXIS, block)
!    do jj = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
!      point(JAXIS) = jj
!      point_mbLp1(JAXIS) = jj - blkLimits(LOW, JAXIS) + 1
!      phys_point(JAXIS) = gr_bhLocCoords(point_mbLp1(JAXIS), JAXIS, block)
!      do ii = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
!        point(IAXIS) = ii
!        point_mbLp1(IAXIS) = ii - blkLimits(LOW, IAXIS) + 1
!        phys_point(IAXIS) = gr_bhLocCoords(point_mbLp1(IAXIS), IAXIS, block)
!        ! insert 1st level blocks into the Priority Queue
!        gr_bhPQSize = 0
!        gr_bhPQSumSquare = 0.0
!        do i = 1, gr_bhTreeNFirstLev
!          tr  = gr_bhTreeFirstLevBlocks(1, i)
!          cpu = gr_bhTreeFirstLevBlocks(2, i)
!          call gr_bhPQInsert(block, tr, cpu, 1, point_mbLp1, phys_point, node)
!        enddo
!   
!        ! walk the tree - replace first elements in queue with their children
!        ! until a desired accuracy is reached
!        do
!   
!          ! pull the first element (interaction) from the queue
!          call gr_bhPQPull(tr, cpu, node_btp, int_type, dr, perr)
!          !print *, "SS after pull: ", sqrt(gr_bhPQSumSquare), perr, tr, cpu, node_btp, gr_bhPQSize, gr_bhTreeMyPE
!
!          select case(int_type)
!          case(GR_TREE_INT_CP) ! CELL-PARENT BLOCK interaction
!            ! insert its children into the queue
!            do l = 1,mchild
!              call gr_bhPQInsert(block, gr_bhTreeChild(1,l,tr,cpu), gr_bhTreeChild(2,l,tr,cpu) &
!                , 1, point_mbLp1, phys_point, node)
!            enddo
!          case(GR_TREE_INT_CN) ! CELL-NODE_OF_BLOCK_TREE interaction
!            if (cpu /= gr_bhTreeMyPE) then
!              level = gr_bhBlockTreePos(GR_TREE_BTP_LEV, node_btp)
!              if ((level+1) > gr_bhLocRecvTreeLevels(tr, cpu)) then
!                print *, "Asking for higher level than present on a given cpu: ", level, gr_bhLocRecvTreeLevels(tr, cpu)
!                print *, "TWB2: ", gr_bhTreeMyPE, block, tr, cpu, sqrt(dr(MDIM+1)) &
!                , gr_bhTreeNodeSize(level+gr_bhTreeLrefine(tr,cpu)) &
!                , gr_bhTreeNodeSize(level+gr_bhTreeLrefine(tr,cpu)) / sqrt(dr(MDIM+1)) &
!                , "|dr:", dr &
!                , "|Cell:", gr_bhLocCoords(point_mbLp1(IAXIS), IAXIS, block), gr_bhLocCoords(point_mbLp1(JAXIS), JAXIS, block), gr_bhLocCoords(point_mbLp1(KAXIS), KAXIS, block)
!                !print *, "TWB3|TA:", gr_bhTreeArray(cpu, tr)%p(pos+GR_TREE_IX-1:pos+GR_TREE_IZ-1)
!                call Driver_abortFlash("TREE_ARRAY ACCESS FAILED!!!")
!              endif
!            endif
!            ! insert its children into the queue
!            do l = GR_TREE_BTP_C1,GR_TREE_BTP_C8
!              call gr_bhPQInsert(block, tr, cpu, gr_bhBlockTreePos(l, node_btp) &
!              & , point_mbLp1, phys_point, node)
!            enddo
!          case default
!            ! DIE  
!            print *, "TreeWalkPQ: all CP and CN interactions eliminated: ", int_type, gr_bhPQSize, sqrt(gr_bhPQSumSquare)
!            ss = 0.0
!            do ipq = 1, gr_bhPQSize
!              print *, "       PQ: ", ipq, gr_bhPriorityQueue(ipq)%tr, gr_bhPriorityQueue(ipq)%perr, ss
!              if (gr_bhPriorityQueue(ipq)%tr /= -1) ss = ss + gr_bhPriorityQueue(ipq)%perr**2
!            enddo
!            print *, "SS before check2: ", sqrt(ss), perr, tr, cpu, node_btp, gr_bhPQSize, gr_bhTreeMyPE
!            call Driver_abortFlash("TreeWalkPQ: all CP and CN interactions eliminated")
!          end select
!
!          !print *, "SS before check1: ", sqrt(gr_bhPQSumSquare), perr, tr, cpu, node_btp, gr_bhPQSize, gr_bhTreeMyPE, gr_bhPriorityQueue(1)%perr
!          !ss = 0.0
!          !do ipq = 1, gr_bhPQSize
!          !  print *, "       PQ: ", ipq, gr_bhPriorityQueue(ipq)%tr, gr_bhPriorityQueue(ipq)%perr, ss
!          !  if (gr_bhPriorityQueue(ipq)%tr /= -1) ss = ss + gr_bhPriorityQueue(ipq)%perr**2
!          !enddo
!          !print *, "SS before check2: ", sqrt(ss), perr, tr, cpu, node_btp, gr_bhPQSize, gr_bhTreeMyPE
!#ifdef GRAVITY
!          if (sqrt(gr_bhPQSumSquare) < grav_bhAccErr) exit
!#endif GRAVITY
!        enddo
!        !print *, "SS final: ", sqrt(gr_bhPQSumSquare), gr_bhTreeMyPE, gr_bhPQSize
!
!        ss = 0.0
!        do ipq = 1, gr_bhPQSize
!          tr       = gr_bhPriorityQueue(ipq)%tr
!          cpu      = gr_bhPriorityQueue(ipq)%cpu
!          node_btp = gr_bhPriorityQueue(ipq)%btp
!          int_type = gr_bhPriorityQueue(ipq)%int_type
!          dr       = gr_bhPriorityQueue(ipq)%dr
!          perr     = gr_bhPriorityQueue(ipq)%perr
!          ss = ss + perr*perr
!          gr_bhPriorityQueue(ipq) = gr_bhPQNull
!          select case(int_type)
!          case(GR_TREE_INT_CP) ! CELL-PARENT BLOCK interaction
!            node(1:GR_TREE_NSIZE) = gr_bhTreeParentTree(1:GR_TREE_NSIZE, tr, cpu)
!            call gr_bhNodeContrib(node, 0, gr_bhTreeLrefine(tr,cpu) &
!            & , dr, block, point, blkLimits, solnData)
!            dcount(3) = dcount(3) + 1
!            !print *, "CP accepted: ", node, "|", gr_bhTreeNodeSize(gr_bhTreeLrefine(tr,cpu)) &
!            !& , "|", dr, "|", block, "|", point, "|", blkLimits
!          case(GR_TREE_INT_CN) ! CELL-NODE_OF_BLOCK_TREE interaction
!            ! leaf block, get the position in gr_bhTreeArray and the node level
!            pos = gr_bhBlockTreePos(GR_TREE_BTP_POS, node_btp)
!            level = gr_bhBlockTreePos(GR_TREE_BTP_LEV, node_btp)
!            node = gr_bhTreeArray(cpu, tr)%p(pos:pos+GR_TREE_NSIZE-1)
!            call gr_bhNodeContrib(node, level, gr_bhTreeLrefine(tr,cpu) &
!            & , dr, block, point, blkLimits, solnData)
!            dcount(2) = dcount(2) + 1
!            !print *, "CN accepted: ", node, "|", gr_bhTreeNodeSize(level+gr_bhTreeLrefine(tr,cpu)) &
!            !& , "|", dr, "|", block, "|", point, "|", blkLimits
!          case(GR_TREE_INT_CC) ! CELL-CELL interaction
!            ! just call BotNodeContrib and thats it
!            pos = gr_bhBlockTreePos(GR_TREE_BTP_POS, node_btp)
!            level = gr_bhBlockTreePos(GR_TREE_BTP_LEV, node_btp)
!            botnode = gr_bhTreeArray(cpu, tr)%p(pos:pos+GR_TREE_BNSIZE-1)
!            call gr_bhBotNodeContrib(botnode, level, gr_bhTreeLrefine(tr,cpu) &
!            & , dr, block, point, blkLimits, solnData)
!            dcount(1) = dcount(1) + 1
!            !print *, "CC accepted: ", botnode, "|", gr_bhTreeNodeSize(level+gr_bhTreeLrefine(tr,cpu)) &
!            !& , "|", dr, "|", block, "|", point, "|", blkLimits
!          case(GR_TREE_INT_SELF)
!            ! add self-contribution of the cell
!            call gr_bhSelfContrib(gr_bhTreeCellSize(gr_bhTreeLrefine(tr,cpu),:) &
!            & , block, point, blkLimits, solnData)
!          case default
!            ! DIE  
!            print *, "TreeWalkPQ: incorrect interaction type: ", int_type
!            call Driver_abortFlash("TreeWalkPQ: incorrect interaction type")
!          end select
!   
!        enddo ! end tree walk
!      enddo
!    enddo
!  enddo
!
!  ! call physics routines to do block finalization
!  call gr_bhFinalizeBlock(block, blkLimits, solnData)
!
!  ! release the block pointer
!  call Grid_releaseBlkPtr(block,solnData, CENTER)
!  gr_bhTreeDcount = gr_bhTreeDcount + dcount
!
!  deallocate(node)
!  deallocate(botnode)
!
!  return
end subroutine gr_bhTreeWalkBlockPQ

