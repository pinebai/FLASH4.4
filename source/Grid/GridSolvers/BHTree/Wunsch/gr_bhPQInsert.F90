!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhPQInsert
!!
!! NAME
!!
!!  gr_bhPQInsert
!!
!!
!! SYNOPSIS
!!
!!  call gr_bhPQInsert(
!!              integer(in)  :: block,
!!              integer(in)  :: tr,
!!              integer(in)  :: cpu,
!!              integer(in)  :: btp,
!!              integer(in)  :: point_mbLp1(MDIM),
!!              real(in)     :: phys_point(MDIM),
!!              real(inout)  :: node(:)
!!             )
!!
!! DESCRIPTION
!!
!!  Inserts node (cpu,tr,btp) into the priority queue.
!!
!! ARGUMENTS
!!
!!  block        : block in which the calculation is done
!!  tr           : block index of the inserted tree node
!!  cpu          : cpu of the inserted tree node
!!  btp          : position of the node of the inserted tree node
!!  point_mbLp1  : point minus blockLimits plus 1 (point-of-calculation)
!!  phys_point   : physical coords of the point-of-calculation
!!  node         : inserted node array
!!
!!
!!***

subroutine gr_bhPQInsert(block, tr, cpu, btp, point_mbLp1, phys_point, node)
  use gr_bhData, ONLY : gr_bhPriorityQueue, gr_bhPQSize, GR_TREE_INT_CC, GR_TREE_INT_CN, &
    GR_TREE_INT_CP, GR_TREE_INT_SELF, gr_bhTreeNodeType, gr_bhTreeParentTree, &
    gr_bhBlockTreePos, gr_bhTreeLevels, gr_bhTreeCellSize, gr_bhTreeMyPE, &
    gr_bhTreeBCen, gr_bhTreeArray, &
    GR_TREE_NSIZE, GR_TREE_BNSIZE, GR_TREE_IX, GR_TREE_IY, GR_TREE_IZ, &
    GR_TREE_BTP_POS ,GR_TREE_BTP_LEV, GR_TREE_BTP_I, GR_TREE_BTP_J, GR_TREE_BTP_K, &
    gr_bhTreeNodeSize, gr_bhTreeBS, gr_bhTreeLrefine, &
    gr_bhTWMaxQueueSize, gr_bhTypePQElement, gr_bhPQSumSquare
  use gr_bhLocalInterface, ONLY : gr_bhPartErr
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  integer, intent(IN)               :: block, tr, cpu, btp
  integer, dimension(MDIM), intent(IN) :: point_mbLp1
  real, dimension(MDIM), intent(IN) :: phys_point
  real, dimension(:), intent(INOUT) :: node
  real, dimension(MDIM+2)           :: dr
  integer :: i, j, k, i_cur, i_par
  integer :: int_type, pos, level
  real :: dist2, perr(2)
  type(gr_bhTypePQElement) :: dummy

  !print *, "PQI called: ", gr_bhTreeMyPE, block, tr, cpu, btp
  !print *, "PQI called2: ", point_mbLp1, phys_point
  ! 1. subtract cell coordinates
  dr(1:MDIM) = -phys_point(:)
  if (gr_bhTreeNodetype(tr, cpu) .ne. 1) then ! 1 = LEAF block
    ! parent block
    int_type = GR_TREE_INT_CP

    ! 2. add position of its mass centre
    dr(IAXIS) = dr(IAXIS) + gr_bhTreeParentTree(GR_TREE_IX, tr, cpu)
    dr(JAXIS) = dr(JAXIS) + gr_bhTreeParentTree(GR_TREE_IY, tr, cpu)
    dr(KAXIS) = dr(KAXIS) + gr_bhTreeParentTree(GR_TREE_IZ, tr, cpu)
    ! fill the node
    node(1:GR_TREE_NSIZE) = gr_bhTreeParentTree(1:GR_TREE_NSIZE, tr, cpu)
  else
    ! leaf block, get the position in gr_bhTreeArray and the node level
    pos = gr_bhBlockTreePos(GR_TREE_BTP_POS, btp)
    level = gr_bhBlockTreePos(GR_TREE_BTP_LEV, btp)

    if (level == gr_bhTreeLevels) then
      ! leaf level of block-tree: cell-cell interaction

      ! get indeces of the interacting cell in the block
      i = gr_bhBlockTreePos(GR_TREE_BTP_I, btp)
      j = gr_bhBlockTreePos(GR_TREE_BTP_J, btp)
      k = gr_bhBlockTreePos(GR_TREE_BTP_K, btp)
            
      ! 2. add position of a leaf node (cell) of the block-tree
      dr(IAXIS) = dr(IAXIS) + gr_bhTreeBCen(IAXIS, tr, cpu) &
      & - 0.5*gr_bhTreeCellSize(gr_bhTreeLrefine(tr, cpu), IAXIS)*(gr_bhTreeBS-2*i+1)
      dr(JAXIS) = dr(JAXIS) + gr_bhTreeBCen(JAXIS, tr, cpu) & 
      & - 0.5*gr_bhTreeCellSize(gr_bhTreeLrefine(tr, cpu), JAXIS)*(gr_bhTreeBS-2*j+1)
      dr(KAXIS) = dr(KAXIS) + gr_bhTreeBCen(KAXIS, tr, cpu) &
      & - 0.5*gr_bhTreeCellSize(gr_bhTreeLrefine(tr, cpu), KAXIS)*(gr_bhTreeBS-2*k+1)

      ! check for self-contribution
      if ((gr_bhTreeMyPE .eq. cpu) .and. (block .eq. tr) &
      & .and. (i .eq. point_mbLp1(IAXIS)) &
      & .and. (j .eq. point_mbLp1(JAXIS)) &
      & .and. (k .eq. point_mbLp1(KAXIS))) then
        int_type = GR_TREE_INT_SELF
      else
        int_type = GR_TREE_INT_CC
      endif

    else
      ! non-leaf nodes of the block tree
      int_type = GR_TREE_INT_CN

      ! 2. add position of the mass centre of the node
      dr(IAXIS) = dr(IAXIS) + gr_bhTreeArray(cpu, tr)%p(pos+GR_TREE_IX-1)
      dr(JAXIS) = dr(JAXIS) + gr_bhTreeArray(cpu, tr)%p(pos+GR_TREE_IY-1)
      dr(KAXIS) = dr(KAXIS) + gr_bhTreeArray(cpu, tr)%p(pos+GR_TREE_IZ-1)
      ! fill the node
      node = gr_bhTreeArray(cpu, tr)%p(pos:pos+GR_TREE_NSIZE-1)
    endif ! test for the leaf level of the block-tree

  endif ! LEAF block or PARENT block
   
  ! 3. correct for periodic boundaries
  call gr_bhPeriodicDr(dr)

  ! calculate square of distance and inverted distance (stored also in dr)
  dist2  = dr(IAXIS)*dr(IAXIS) + dr(JAXIS)*dr(JAXIS) + dr(KAXIS)*dr(KAXIS)
  dr(MDIM+1) = dist2
  dr(MDIM+2) = sqrt(1.0 / (dist2 + 1D-99))

  ! calculate maximum partial error of this interaction
  if ((int_type == GR_TREE_INT_CP) .or. (int_type == GR_TREE_INT_CN)) then
    call gr_bhPartErr(node, gr_bhTreeNodeSize(gr_bhTreeLrefine(tr,cpu)), dr, perr)
  else
    perr = 0.0
  endif
  gr_bhPQSumSquare = gr_bhPQSumSquare + perr(1)*perr(1)
  !print *, "PQI: added to SS: ", gr_bhTreeMyPE, perr(1)**2

  ! write data into the Priority Queue array
  gr_bhPQSize = gr_bhPQSize + 1
  i_cur = gr_bhPQSize ! index of a current element
  gr_bhPriorityQueue(i_cur)%tr       = tr
  gr_bhPriorityQueue(i_cur)%cpu      = cpu
  gr_bhPriorityQueue(i_cur)%btp      = btp
  gr_bhPriorityQueue(i_cur)%int_type = int_type
  gr_bhPriorityQueue(i_cur)%dr       = dr
  gr_bhPriorityQueue(i_cur)%perr     = perr(1)

  ! binary heap insert - ensures that the element with the largest perr 
  ! will be the first one in the queue
  do
    i_par = i_cur / 2
    if (i_par > 0) then ! fortran does not guarantee short-circuit evaluation :(
      if (gr_bhPriorityQueue(i_par)%perr < gr_bhPriorityQueue(i_cur)%perr) then
        dummy = gr_bhPriorityQueue(i_par)
        gr_bhPriorityQueue(i_par) = gr_bhPriorityQueue(i_cur)
        gr_bhPriorityQueue(i_cur) = dummy
        i_cur = i_par
      else
        exit
      endif
    else
      exit
    endif
  enddo

  !do i = 2,gr_bhTWMaxQueueSize
  !  if (gr_bhPriorityQueue(i)%perr > gr_bhPriorityQueue(1)%perr) then
  !    ! DIE
  !    print *, "PQInsert: PQ is not sorted after insert"
  !    call Driver_abortFlash("PQInsert: PQ is not sorted after insert")
  !  endif
  !enddo


  !print *, "PQI end: ", gr_bhPriorityQueue(gr_bhPQSize)%int_type, gr_bhPriorityQueue(gr_bhPQSize)%dr, gr_bhPriorityQueue(gr_bhPQSize)%perr
  return
end subroutine gr_bhPQInsert


