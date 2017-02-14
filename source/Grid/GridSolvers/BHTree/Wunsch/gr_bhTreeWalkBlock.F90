!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhTreeWalkBlock
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
!!   Traverses the tree in a specific block. In this, "standard tree walk",
!!   the MAC for cell-cell and cell-tree_node interactions is evaluated for 
!!   each cell, however, the MAC for cell-parent_tree interaction is evaluated
!!   for whole block just once.
!!
!! ARGUMENTS
!!
!!  block    - ID of a block where the contribution is calculated
!!
!!***



subroutine gr_bhTreeWalkBlock(block)

  use Grid_interface, ONLY : Grid_getLocalNumBlks, &
    Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_bhLocalInterface, ONLY : gr_bhLeafContrib, gr_bhParentContrib &
    , gr_bhStartBlock, gr_bhFinalizeBlock, gr_bhMAC, gr_bhPeriodicDr
  use gr_bhData, ONLY : gr_bhTreeLevels, &
    gr_bhTreeBS, gr_bhLocCoords, gr_bhTreeMyPE, &
    gr_bhTreeLrefine, gr_bhTreeNodeSize2, gr_bhPhysMACTW_step, &
    gr_bhTreeDcount, gr_bhTreeLimangle2i, &
    gr_bhTreeNFirstLev, gr_bhTreeChild, gr_bhTreeNodetype, &
    gr_bhTreeFirstLevBlocks, gr_bhTreeParentTree, &
    GR_TREE_IX, GR_TREE_IY, GR_TREE_IZ, GR_TREE_NSIZE, &
    gr_bhTreeCellSize, gr_bhTreeBCen
  use tree, ONLY : nodetype, child, mchild

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  integer, intent(IN) :: block
  integer       :: i, cpu, tr, sp, ii, istat
  real          :: dcount(1:3), cellcnt, nodecnt
  real, POINTER, DIMENSION(:,:,:,:) :: solnData
  integer :: stack(2, MAXBLOCKS)
  real    :: dist2
  integer, dimension(2,MDIM)   :: blkLimits,blkLimitsGC
  real, dimension(MDIM+2) :: dr ! two additional elements dr^2 and 1/|dr|
  real, dimension(MDIM) :: drDiff, drGC
  integer, dimension(MDIM) :: point
  real, allocatable  :: node(:)
  logical       :: node_accepted, gcell = .false.

  ! allocate the node array
  allocate(node(GR_TREE_NSIZE), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate node in gr_bhTreeWalkBlock.F90")

  ! get information about the block
  call Grid_getBlkIndexLimits(block,blkLimits,blkLimitsGC) ! indexes necessary to write values into solnData
  call Grid_getBlkPtr(block,solnData,CENTER) ! pointer to the density and gp field

  ! call physics routines to do block initialization
  call gr_bhStartBlock(block, blkLimits, solnData)
  
  ! reset distance counter
  dcount = 0.0

  ! put 1st level blocks on stack
  stack(1:2,1:gr_bhTreeNFirstLev) = gr_bhTreeFirstLevBlocks(1:2,1:gr_bhTreeNFirstLev)
  sp = gr_bhTreeNFirstLev

  ! walk the amr tree
  do

    ! pop a block from the stack
    tr = stack(1,sp)
    cpu  = stack(2,sp)
    sp = sp - 1

    if (gr_bhTreeNodetype(tr, cpu) .eq. 1) then ! 1 = LEAF block
      ! contribution of a LEAF block
      call gr_bhLeafContrib(block, tr, cpu, solnData, blkLimits, cellcnt, nodecnt)
      dcount(1) = dcount(1) + cellcnt
      dcount(2) = dcount(2) + nodecnt
    else
      ! contribution of a parent block
 
      ! determine Bmin vector: distance between the closest cell of this block 
      ! and the closest cell of the parent block (tr, cpu)
      ! 1. distance between this block and the parent block centres
      dr(IAXIS) = gr_bhTreeParentTree(GR_TREE_IX, tr, cpu) &
      &         - gr_bhLocCoords(gr_bhTreeBS+1, IAXIS, block)
      dr(JAXIS) = gr_bhTreeParentTree(GR_TREE_IY, tr, cpu) &
      &         - gr_bhLocCoords(gr_bhTreeBS+1, JAXIS, block)
      dr(KAXIS) = gr_bhTreeParentTree(GR_TREE_IZ, tr, cpu) &
      &         - gr_bhLocCoords(gr_bhTreeBS+1, KAXIS, block)
      ! 2. Correct dr for periodic boundaries
      call gr_bhPeriodicDr(dr)
      ! 3. subtract this block and parrent block sizes, make sure the vector
      ! component does not change the sign
      if (dr(IAXIS) < 0) then
        dr(IAXIS) = min(0.0, dr(IAXIS) + 0.5 * gr_bhTreeBS * ( &
        & gr_bhTreeCellSize(gr_bhTreeLrefine(block, gr_bhTreeMyPE), IAXIS) + &
        & gr_bhTreeCellSize(gr_bhTreeLrefine(tr, cpu), IAXIS)))
      else
        dr(IAXIS) = max(0.0, dr(IAXIS) - 0.5 * gr_bhTreeBS * ( &
        & gr_bhTreeCellSize(gr_bhTreeLrefine(block, gr_bhTreeMyPE), IAXIS) + &
        & gr_bhTreeCellSize(gr_bhTreeLrefine(tr, cpu), IAXIS)))
      endif
      if (dr(JAXIS) < 0) then
        dr(JAXIS) = min(0.0, dr(JAXIS) + 0.5 * gr_bhTreeBS * ( &
        & gr_bhTreeCellSize(gr_bhTreeLrefine(block, gr_bhTreeMyPE), JAXIS) + &
        & gr_bhTreeCellSize(gr_bhTreeLrefine(tr, cpu), JAXIS)))
      else
        dr(JAXIS) = max(0.0, dr(JAXIS) - 0.5 * gr_bhTreeBS * ( &
        & gr_bhTreeCellSize(gr_bhTreeLrefine(block, gr_bhTreeMyPE), JAXIS) + &
        & gr_bhTreeCellSize(gr_bhTreeLrefine(tr, cpu), JAXIS)))
      endif
      if (dr(KAXIS) < 0) then
        dr(KAXIS) = min(0.0, dr(KAXIS) + 0.5 * gr_bhTreeBS * ( &
        & gr_bhTreeCellSize(gr_bhTreeLrefine(block, gr_bhTreeMyPE), KAXIS) + &
        & gr_bhTreeCellSize(gr_bhTreeLrefine(tr, cpu), KAXIS)))
      else
        dr(KAXIS) = max(0.0, dr(KAXIS) - 0.5 * gr_bhTreeBS * ( &
        & gr_bhTreeCellSize(gr_bhTreeLrefine(block, gr_bhTreeMyPE), KAXIS) + &
        & gr_bhTreeCellSize(gr_bhTreeLrefine(tr, cpu), KAXIS)))
      endif

      ! construct the node array
      node(1:GR_TREE_NSIZE) = gr_bhTreeParentTree(1:GR_TREE_NSIZE, tr, cpu)

      ! position vector to the node mass centre
      dist2  = dr(IAXIS)*dr(IAXIS) + dr(JAXIS)*dr(JAXIS) + dr(KAXIS)*dr(KAXIS)
      dr(MDIM+1) = dist2
      dr(MDIM+2) = sqrt(1.0 / (dist2 + 1D-99))

      ! position vector to the node geometric centre
      ! drDiff is a vector from the node mass centre to its geometrical centre
      drDiff(:) = gr_bhTreeBCen(:, tr, cpu) - node(GR_TREE_IX:GR_TREE_IZ)
      !drDiff(IAXIS) = gr_bhTreeBCen(IAXIS, tr, cpu) - node(GR_TREE_IX)
      !drDiff(JAXIS) = gr_bhTreeBCen(JAXIS, tr, cpu) - node(GR_TREE_IY)
      !drDiff(KAXIS) = gr_bhTreeBCen(KAXIS, tr, cpu) - node(GR_TREE_IZ)
      drGC = dr(1:MDIM) + drDiff

      ! the node acceptance criterion (MAC)
      point(:) = -1
      node_accepted = gr_bhMAC(gr_bhPhysMACTW_step, node &
      & , gr_bhTreeNodeSize2(gr_bhTreeLrefine(tr,cpu)) &
      & , drGC, dr, block, point, blkLimits, solnData)

      if (node_accepted) then 
        ! node accepted, call contribution subroutines
        call gr_bhParentContrib(block, tr, cpu, solnData, blkLimits)
        ! increase the distance counter by a number of cells in this block
        dcount(3) = dcount(3) + gr_bhTreeBS**3 
      else
        ! node not accepted, put all its children on the stack
        do ii = 1,mchild
          sp = sp + 1
          stack(1,sp) = gr_bhTreeChild(1,ii,tr,cpu)
          stack(2,sp) = gr_bhTreeChild(2,ii,tr,cpu)
        enddo

      endif ! add or use children?
      
    endif ! LEAF block


    ! stack is empty - exiting
    if (sp == 0) exit
  enddo

  ! call physics routines to do block finalization
  call gr_bhFinalizeBlock(block, blkLimits, solnData)

  ! release the block pointer
  call Grid_releaseBlkPtr(block,solnData, CENTER)
  gr_bhTreeDcount = gr_bhTreeDcount + dcount

  deallocate(node)

  return
end subroutine gr_bhTreeWalkBlock

